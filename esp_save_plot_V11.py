#!/usr/bin/env python3
"""
ESP32 12-Lead ECG Real-Time WiFi Plotter  (v11)
=============================================================
Connects to ESP32 AP on TCP port 3333, receives binary-framed
ECG data (8 channels @ 500 SPS, RDATA mode), computes derived
leads (III, aVR, aVL, aVF) for full 12-lead display.

Features:
  • ADS1298 status byte validation (0xCx) – bad samples skipped
  • Full 12-lead ECG display (I, II, III, aVR, aVL, aVF, V1–V6)
  • OpenGL-accelerated rendering (~60 FPS)
  • CSV recording with full 500 SPS resolution
  • Remote PGA gain control via TCP command port
  • Real-time HR estimation (threshold-based R-peak detection)
  • Double-click to maximize/restore individual channel

Binary protocol (must match main.cpp v11):
    Frame = [0xA5][0x5A][count][sample0]...[sampleN-1]
    Sample = [4B timestamp_us LE][27B raw ADS1298 SPI data]
    27B = 3B status + 8×3B channels (24-bit signed, big-endian)

Usage:
    1. Connect your PC WiFi to "ESP32_ECG" (password: 12345678)
    2. Run: python esp_save_plot_V11.py
    3. The plotter auto-connects and starts streaming

Dependencies:
    pip install pyqtgraph PyQt5 numpy
"""

import csv
import os
import socket
import struct
import threading
import time
import sys
from datetime import datetime

import numpy as np
import pyqtgraph as pg
from PyQt5 import QtWidgets, QtCore, QtGui

# =============================================================================
# CONFIGURATION
# =============================================================================

ESP_HOST        = "192.168.4.1"     # ESP32 AP default IP
TCP_PORT        = 3333              # Data port
TCP_PORT_CMD    = 3334              # Command port (PC → ESP)
SAMPLE_RATE     = 500.0             # Hz (ADS1298 @ 500 SPS)
NUM_RAW_CH      = 8                 # 8 ADS1298 raw channels
NUM_DERIVED     = 4                 # III, aVR, aVL, aVF
NUM_CHANNELS    = NUM_RAW_CH + NUM_DERIVED   # 12 total

# Ring buffer / display
WINDOW_SEC      = 10.0               # Seconds of data visible on screen
WINDOW_SAMPLES  = int(SAMPLE_RATE * WINDOW_SEC)

# Plot refresh interval (ms) — menší = plynulejší, 12 ms ≈ 83 FPS
PLOT_INTERVAL_MS = 12
# Plynulý posun osy X: 0 = žádné, 1 = okamžitý skok; 0.2–0.35 = plynulé sledování
SCROLL_SMOOTH_ALPHA = 0.28

# # --- Y_SPAN: zakomentováno kvůli manuálnímu škálování kolečkem myši ---
# # Minimální rozsah osy Y (mV) pro jeden svod; skutečný rozsah = max(rozsah dat, Y_SPAN_MIN_MV)
# # → každý kanál se přizpůsobí amplitudě svého signálu, data neklipují
# Y_SPAN_MIN_MV = 0.5
# # Volitelně: pevný rozsah pro každý svod (mV), None = auto z dat
# # Příklad: [2.0, 2.0, 1.5, 1.5, 1.5, 1.5, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0] pro I, II, III, aVR, aVL, aVF, V1–V6
# Y_SPAN_PER_CHANNEL_MV = None  # nebo seznam 12 čísel

# TCP receive buffer size
TCP_RECV_SIZE   = 16384

# ---- Binary protocol constants (must match ESP32 firmware v11) ----
SYNC_MARKER     = bytes([0xA5, 0x5A])
RAW_SPI_BYTES   = 27          # 3B status + 8×3B channels
SAMPLE_SIZE     = 4 + RAW_SPI_BYTES   # 31 bytes: 4B ts + 27B raw
FRAME_HEADER    = 3           # 2B sync + 1B count
SAMPLES_PER_FRAME = 10
FRAME_EXPECTED  = FRAME_HEADER + SAMPLES_PER_FRAME * SAMPLE_SIZE  # 623

# ---- ADS1298 → Voltage conversion ----
# V_input = ADC_code × VREF / (Gain × 2^23)
# With VREF=2.4 V, Gain=6:  1 LSB = 4.76837158e-8 V = 0.0476837158 µV
ADS_VREF        = 2.4              # Internal reference voltage [V]
ADS_GAIN        = 6                # PGA gain (CHnSET[6:4] = 000)
ADS_LSB_MV      = (ADS_VREF / (ADS_GAIN * (2**23))) * 1e3   # mV per LSB

# Valid ADS1298 PGA gain values
VALID_GAINS     = [1, 2, 3, 4, 6, 8, 12]

def compute_lsb_mv(gain: int) -> float:
    """Compute mV-per-LSB for a given PGA gain."""
    return (ADS_VREF / (gain * (2**23))) * 1e6 / 1000.0

# 12-lead display order:  I, II, III, aVR, aVL, aVF, V1–V6
# Index mapping from raw ADS1298 data (0–7) to display channels (0–11):
#   raw CH1=V1, CH2=Lead I, CH3=Lead II, CH4=V2, CH5=V3, CH6=V4, CH7=V5, CH8=V6
#   derived: III = II−I,  aVR = −(I+II)/2,  aVL = I−II/2,  aVF = II−I/2
CHANNEL_LABELS = [
    "I  (LA−RA)",  "II (LL−RA)",  "III",  "aVR",  "aVL",  "aVF",
    "V1",  "V2",  "V3",  "V4",  "V5",  "V6",
]

CHANNEL_COLORS = [
    (  0, 200,   0),   # I    – green
    (  0, 120, 255),   # II   – blue
    (  0, 180, 180),   # III  – teal
    (200, 100,   0),   # aVR  – brown
    (180,   0, 255),   # aVL  – purple
    (255, 165,   0),   # aVF  – orange
    (255,  80,  80),   # V1   – red
    (255, 200,  40),   # V2   – yellow
    (100, 255, 100),   # V3   – lime
    (  0, 210, 210),   # V4   – cyan
    (255, 100, 200),   # V5   – pink
    (200, 200, 200),   # V6   – white
]


# Directory where this script lives – CSV files are saved here by default
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
LOGS_DIR   = os.path.join(SCRIPT_DIR, "logs")

# CSV column headers (matches existing log format)
CSV_HEADER = [
    "timestamp_us",
    "I_mV", "II_mV", "III_mV", "aVR_mV", "aVL_mV", "aVF_mV",
    "V1_mV", "V2_mV", "V3_mV", "V4_mV", "V5_mV", "V6_mV",
    "gain",
]


# =============================================================================
# CSV RECORDER  (thread-safe, called from receiver thread)
# =============================================================================

class CsvRecorder:
    """Thread-safe CSV writer.  Called from TcpReceiver thread so every
    decoded sample at full 500 SPS is written – completely independent of
    the display ring-buffer which is windowed / down-sampled."""

    def __init__(self):
        self._lock = threading.Lock()
        self._file = None
        self._writer = None
        self._sample_count = 0
        self._filename = ""
        self._recording = False
        self._start_time = 0.0

    # -- public API (called from GUI thread) --------------------------------
    @property
    def recording(self) -> bool:
        with self._lock:
            return self._recording

    @property
    def sample_count(self) -> int:
        with self._lock:
            return self._sample_count

    @property
    def filename(self) -> str:
        with self._lock:
            return self._filename

    @property
    def elapsed(self) -> float:
        with self._lock:
            if not self._recording:
                return 0.0
            return time.time() - self._start_time

    def start(self, filepath: str | None = None):
        """Open a new CSV file and start recording."""
        with self._lock:
            if self._recording:
                return
            if filepath is None:
                os.makedirs(LOGS_DIR, exist_ok=True)
                stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                filepath = os.path.join(LOGS_DIR, f"ecg_data_{stamp}.csv")
            self._filename = filepath
            self._file = open(filepath, "w", newline="", buffering=8192)
            self._writer = csv.writer(self._file)
            self._writer.writerow(CSV_HEADER)
            self._sample_count = 0
            self._start_time = time.time()
            self._recording = True

    def stop(self):
        """Flush and close the CSV file."""
        with self._lock:
            if not self._recording:
                return
            self._recording = False
            try:
                self._file.flush()
                self._file.close()
            except Exception:
                pass
            self._file = None
            self._writer = None

    # -- called from receiver thread ----------------------------------------
    def write_samples(self, timestamps, channels, gain: int):
        """Write a batch of samples.  *timestamps*: list[int] (µs),
        *channels*: list[list[float]] (12 values in mV each)."""
        with self._lock:
            if not self._recording or self._writer is None:
                return
            for ts, ch in zip(timestamps, channels):
                row = [ts] + [f"{v:.6f}" for v in ch] + [gain]
                self._writer.writerow(row)
            self._sample_count += len(timestamps)
            # Flush every ~500 samples (≈1 s) to keep file up-to-date
            if self._sample_count % 500 < len(timestamps):
                try:
                    self._file.flush()
                except Exception:
                    pass


# =============================================================================
# RING BUFFER  (thread-safe, numpy-backed)
# =============================================================================

class RingBuffer:
    """Fixed-capacity ring buffer for multi-channel float data + timestamps."""

    def __init__(self, capacity: int, num_channels: int):
        self.capacity = capacity
        self.num_ch = num_channels
        self.lock = threading.Lock()

        # shape: (num_channels, capacity)
        self.data = np.zeros((num_channels, capacity), dtype=np.float32)
        # timestamps in microseconds (uint64 to avoid float precision loss)
        self.timestamps = np.zeros(capacity, dtype=np.uint64)
        self.head = 0           # next write position
        self.count = 0          # valid samples (≤ capacity)

    def push(self, timestamps: np.ndarray, samples: np.ndarray):
        """Push N samples.
        timestamps shape: (N,) — uint64 microseconds
        samples shape: (N, num_channels)."""
        n = samples.shape[0]
        if n == 0:
            return
        with self.lock:
            cap = self.capacity
            # Keep only newest if more than capacity
            if n > cap:
                timestamps = timestamps[-cap:]
                samples = samples[-cap:]
                n = cap

            h = self.head
            first = min(n, cap - h)
            second = n - first

            # Bulk-copy transposed to (num_ch, cap)
            self.data[:, h:h + first] = samples[:first].T
            self.timestamps[h:h + first] = timestamps[:first]
            if second > 0:
                self.data[:, 0:second] = samples[first:].T
                self.timestamps[0:second] = timestamps[first:]

            self.head = (h + n) % cap
            self.count = min(self.count + n, cap)

    def snapshot(self):
        """Return time-ordered copy: (time_array_sec, data_array, timestamps_us).
        time_array_sec is relative: newest = 0, oldest negative.
        data shape: (num_channels, count)."""
        with self.lock:
            n = self.count
            h = self.head
            if n == 0:
                return np.array([]), np.empty((self.num_ch, 0)), np.array([], dtype=np.uint64)

            if n < self.capacity:
                d = self.data[:, :n].copy()
                ts = self.timestamps[:n].copy()
            else:
                # Reorder wrapped buffer to chronological order
                d = np.empty((self.num_ch, self.capacity), dtype=np.float32)
                ts = np.empty(self.capacity, dtype=np.uint64)
                tail = self.capacity - h
                d[:, :tail] = self.data[:, h:]
                d[:, tail:] = self.data[:, :h]
                ts[:tail] = self.timestamps[h:]
                ts[tail:] = self.timestamps[:h]

            # Use ESP32 timestamps to build time axis (seconds since ESP32 boot)
            t = ts.astype(np.float64) / 1_000_000.0
            t = t.astype(np.float32)

            return t, d, ts

    def clear(self):
        with self.lock:
            self.data.fill(0)
            self.timestamps.fill(0)
            self.head = 0
            self.count = 0


# =============================================================================
# TCP RECEIVER THREAD
# =============================================================================

class TcpReceiver(threading.Thread):
    """Background thread: connects to ESP32, receives binary frames,
    decodes raw ADS1298 data, computes derived leads, validates status
    bytes, and pushes 12-lead samples into ring buffer."""

    def __init__(self, ring: RingBuffer, recorder: CsvRecorder | None = None):
        super().__init__(daemon=True)
        self.ring = ring
        self.recorder = recorder
        self.running = True
        self.connected = False
        self.packets_received = 0
        self.errors = 0
        self.bad_status_count = 0
        self._status = "Odpojeno"
        self._lock = threading.Lock()
        self.last_latency_us = 0
        self.frames_received = 0
        self.current_gain = ADS_GAIN
        self.lsb_mv = ADS_LSB_MV

    @property
    def status(self):
        with self._lock:
            return self._status

    @status.setter
    def status(self, val):
        with self._lock:
            self._status = val

    def stop(self):
        self.running = False

    def _decode_frames(self, buf: bytes):
        """Parse binary frames from byte buffer.
        Returns (remaining_bytes, timestamps_list, channels_list).
        Each decoded sample yields one timestamp (uint32) and 12 float channels
        (I, II, III, aVR, aVL, aVF, V1–V6) with status validation."""
        results_ts = []
        results_ch = []
        frames_in_batch = 0
        pos = 0
        buf_len = len(buf)

        while pos <= buf_len - FRAME_HEADER:
            # Search for sync marker
            sync_pos = buf.find(SYNC_MARKER, pos)
            if sync_pos == -1:
                pos = max(pos, buf_len - 1)
                break

            pos = sync_pos
            if pos + FRAME_HEADER > buf_len:
                break

            count = buf[pos + 2]
            if count == 0 or count > SAMPLES_PER_FRAME:
                pos += 2
                continue

            frame_end = pos + FRAME_HEADER + count * SAMPLE_SIZE
            if frame_end > buf_len:
                break  # incomplete frame, wait for more data

            # Decode all samples in this frame
            frames_in_batch += 1
            for s in range(count):
                off = pos + FRAME_HEADER + s * SAMPLE_SIZE

                # --- Status byte validation (ADS1298: top nibble = 0xC) ---
                status_byte = buf[off + 4]  # first byte after 4B timestamp
                if (status_byte & 0xF0) != 0xC0:
                    self.bad_status_count += 1
                    continue  # skip this sample entirely

                # Timestamp: 4 bytes little-endian uint32
                ts_us = struct.unpack_from('<I', buf, off)[0]

                # Raw ADS1298: 3B status + 8×3B channel data
                raw_ch = []
                for c in range(NUM_RAW_CH):
                    idx = off + 4 + 3 + c * 3   # skip ts(4) + status(3)
                    val = (buf[idx] << 16) | (buf[idx + 1] << 8) | buf[idx + 2]
                    if val & 0x800000:
                        val -= 0x1000000
                    raw_ch.append(val)

                # --- Compute derived 12-lead channels ---
                # raw_ch indices: 0=V1, 1=LeadI(CH2), 2=LeadII(CH3),
                #                 3=V2, 4=V3, 5=V4, 6=V5, 7=V6
                lead_I  = raw_ch[1]
                lead_II = raw_ch[2]
                lead_III = lead_II - lead_I
                aVR = -(lead_I + lead_II) / 2.0
                aVL = lead_I - lead_II / 2.0
                aVF = lead_II - lead_I / 2.0

                lsb = self.lsb_mv  # use current gain-dependent LSB

                # 12-lead order: I, II, III, aVR, aVL, aVF, V1, V2, V3, V4, V5, V6
                # Convert from raw ADC counts → millivolts
                ch12 = [
                    float(lead_I)   * lsb,
                    float(lead_II)  * lsb,
                    float(lead_III) * lsb,
                    aVR             * lsb,
                    aVL             * lsb,
                    aVF             * lsb,
                    float(raw_ch[0]) * lsb,   # V1
                    float(raw_ch[3]) * lsb,   # V2
                    float(raw_ch[4]) * lsb,   # V3
                    float(raw_ch[5]) * lsb,   # V4
                    float(raw_ch[6]) * lsb,   # V5
                    float(raw_ch[7]) * lsb,   # V6
                ]

                results_ts.append(ts_us)
                results_ch.append(ch12)

            pos = frame_end

        return buf[pos:], results_ts, results_ch, frames_in_batch

    def run(self):
        recv_buf = b""
        while self.running:
            # --- connect ---
            try:
                self.status = f"Připojování k {ESP_HOST}:{TCP_PORT} …"
                sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                sock.settimeout(5.0)
                sock.connect((ESP_HOST, TCP_PORT))
                sock.setsockopt(socket.IPPROTO_TCP, socket.TCP_NODELAY, 1)
                sock.settimeout(2.0)
                self.connected = True
                self.status = "Připojeno ✓"
                recv_buf = b""
                self._esp_time_offset = None
            except Exception as e:
                self.status = f"Chyba spojení: {e}"
                self.connected = False
                time.sleep(2)
                continue

            # --- receive loop ---
            while self.running and self.connected:
                try:
                    raw = sock.recv(TCP_RECV_SIZE)
                    if not raw:
                        self.connected = False
                        self.status = "Spojení uzavřeno"
                        break

                    recv_time_us = int(time.time() * 1_000_000)
                    recv_buf += raw

                    # Decode as many complete frames as possible
                    recv_buf, batch_ts, batch_ch, n_frames = self._decode_frames(recv_buf)

                    if batch_ts:
                        ts_arr = np.array(batch_ts, dtype=np.uint64)
                        ch_arr = np.array(batch_ch, dtype=np.float32)
                        self.ring.push(ts_arr, ch_arr)
                        self.packets_received += len(batch_ts)
                        self.frames_received += n_frames

                        # Write to CSV recorder (full rate, no downsampling)
                        if self.recorder is not None:
                            self.recorder.write_samples(
                                batch_ts, batch_ch, self.current_gain
                            )

                        # Estimate one-way latency from last sample in batch
                        if self._esp_time_offset is None:
                            self._esp_time_offset = recv_time_us - int(ts_arr[-1])
                        self.last_latency_us = recv_time_us - (int(ts_arr[-1]) + self._esp_time_offset)

                    # Safety: if recv_buf grows too large without valid frames, trim it
                    if len(recv_buf) > TCP_RECV_SIZE * 4:
                        recv_buf = recv_buf[-FRAME_EXPECTED:]
                        self.errors += 1

                except socket.timeout:
                    continue
                except Exception as e:
                    self.status = f"Chyba příjmu: {e}"
                    self.connected = False
                    self.errors += 1
                    break

            # --- cleanup ---
            try:
                sock.close()
            except:
                pass
            if self.running:
                self.status = "Odpojeno – opakuji za 2 s …"
                time.sleep(2)


# =============================================================================
# REAL-TIME PLOT WINDOW  (PyQtGraph)
# =============================================================================

class ECGPlotWindow(QtWidgets.QMainWindow):
    """Main window with stacked ECG channel subplots and a status bar."""

    def __init__(self, ring: RingBuffer, receiver: TcpReceiver, recorder: CsvRecorder):
        super().__init__()
        self.ring = ring
        self.receiver = receiver
        self.recorder = recorder

        self.setWindowTitle("ESP32 12-Lead ECG  —  Real-Time WiFi Plotter (v11)")
        self.resize(1400, 1100)

        # --- central widget ---
        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        layout = QtWidgets.QVBoxLayout(central)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(2)

        # --- pyqtgraph layout ---
        pg.setConfigOptions(antialias=False, useOpenGL=True)
        self.pw = pg.GraphicsLayoutWidget()
        self.pw.setBackground("k")
        layout.addWidget(self.pw, stretch=1)

        # Create one subplot per channel (12 leads), stacked vertically
        self.plots = []
        self.curves = []
        for i in range(NUM_CHANNELS):
            if i > 0:
                self.pw.nextRow()
            p = self.pw.addPlot()
            p.setLabel("left", CHANNEL_LABELS[i], units="mV")
            # Stejná šířka prostoru pro popisky Y u všech grafů → osy vizuálně srovnané
            p.getAxis("left").setStyle(tickTextWidth=52, autoExpandTextSpace=False)
            p.showGrid(x=True, y=True, alpha=0.25)
            p.setMouseEnabled(x=False, y=True)
            p.enableAutoRange(axis="y", enable=True)
            p.setClipToView(True)
            p.setDownsampling(mode='peak')
            # Compact vertical height
            p.setMinimumHeight(40)
            p.setMaximumHeight(120)
            # Hide x-axis labels except last subplot
            if i < NUM_CHANNELS - 1:
                p.getAxis("bottom").setStyle(showValues=False)
            else:
                p.setLabel("bottom", "Čas", units="s")

            # Link X axes for synchronised zoom/pan
            if i > 0:
                p.setXLink(self.plots[0])

            pen = pg.mkPen(color=CHANNEL_COLORS[i], width=1.2)
            curve = p.plot(pen=pen)
            self.plots.append(p)
            self.curves.append(curve)

        # --- double-click to maximize/restore ---
        self._maximized_index = None   # None = all visible, int = maximized channel
        self.pw.scene().sigMouseClicked.connect(self._on_scene_click)

        # --- bottom controls bar ---
        controls_layout = QtWidgets.QHBoxLayout()
        controls_layout.setContentsMargins(4, 0, 4, 0)

        # Gain selector
        gain_label = QtWidgets.QLabel("GAIN:")
        gain_label.setStyleSheet("color: #ccc; font-weight: bold; font-size: 13px;")
        controls_layout.addWidget(gain_label)

        self.gain_combo = QtWidgets.QComboBox()
        self.gain_combo.setFixedWidth(80)
        for g in VALID_GAINS:
            self.gain_combo.addItem(f"x{g}", g)
        # Set default to current gain (6)
        default_idx = VALID_GAINS.index(ADS_GAIN)
        self.gain_combo.setCurrentIndex(default_idx)
        self.gain_combo.setStyleSheet(
            "QComboBox { background: #333; color: #0f0; font-size: 13px; "
            "border: 1px solid #555; padding: 2px 6px; }"
            "QComboBox::drop-down { border: none; }"
            "QComboBox QAbstractItemView { background: #333; color: #0f0; "
            "selection-background-color: #555; }"
        )
        self.gain_combo.currentIndexChanged.connect(self._on_gain_changed)
        controls_layout.addWidget(self.gain_combo)

        self.gain_status_label = QtWidgets.QLabel("")
        self.gain_status_label.setStyleSheet(
            "color: #888; font-family: Consolas, monospace; font-size: 12px; padding-left: 8px;"
        )
        controls_layout.addWidget(self.gain_status_label)

        # --- Spacer ---
        controls_layout.addSpacing(24)

        # --- REC button ---
        self.rec_btn = QtWidgets.QPushButton("⏺  REC")
        self.rec_btn.setFixedWidth(100)
        self.rec_btn.setCheckable(True)
        self.rec_btn.setStyleSheet(
            "QPushButton { background: #333; color: #ccc; font-size: 13px; "
            "font-weight: bold; border: 1px solid #555; padding: 2px 8px; }"
            "QPushButton:checked { background: #800; color: #f44; border: 1px solid #f44; }"
        )
        self.rec_btn.clicked.connect(self._on_rec_toggled)
        controls_layout.addWidget(self.rec_btn)

        self.rec_status_label = QtWidgets.QLabel("")
        self.rec_status_label.setStyleSheet(
            "color: #888; font-family: Consolas, monospace; font-size: 12px; padding-left: 8px;"
        )
        controls_layout.addWidget(self.rec_status_label)

        # --- Spacer ---
        controls_layout.addSpacing(24)

        # --- Pozastavit / Spustit vykreslování ---
        self._plot_paused = False
        self.pause_btn = QtWidgets.QPushButton("⏸ Pozastavit vykreslování")
        self.pause_btn.setFixedWidth(180)
        self.pause_btn.setCheckable(True)
        self.pause_btn.setStyleSheet(
            "QPushButton { background: #333; color: #ccc; font-size: 13px; "
            "font-weight: bold; border: 1px solid #555; padding: 2px 8px; }"
            "QPushButton:checked { background: #084; color: #4f4; border: 1px solid #4f4; }"
        )
        self.pause_btn.clicked.connect(self._on_pause_toggled)
        controls_layout.addWidget(self.pause_btn)

        controls_layout.addStretch()
        layout.addLayout(controls_layout)

        # --- bottom status bar ---
        self.status_label = QtWidgets.QLabel("Spouštím …")
        self.status_label.setStyleSheet(
            "color: #aaa; font-family: Consolas, monospace; font-size: 12px; padding: 2px;"
        )
        layout.addWidget(self.status_label)

        # --- plynulý posun osy X (aktuální rozsah viewu) ---
        self._x_view_min = None
        self._x_view_max = None

        # --- refresh timer ---
        self.timer = QtCore.QTimer()
        self.timer.timeout.connect(self._update)
        self.timer.start(PLOT_INTERVAL_MS)

        self._last_stat_time = time.time()
        self._last_pkt_count = 0
        self._last_hr_time = 0.0
        self._last_hr_str = "HR: —"

        # --- command socket (for gain etc.) ---
        self._cmd_sock = None
        self._cmd_lock = threading.Lock()

        # --- recording status update timer ---
        self._rec_timer = QtCore.QTimer()
        self._rec_timer.timeout.connect(self._update_rec_status)
        self._rec_timer.start(500)

    # ------------------------------------------------------------------ #
    def _connect_cmd(self):
        """Ensure command TCP socket is connected."""
        with self._cmd_lock:
            if self._cmd_sock is not None:
                return True
            try:
                s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                s.settimeout(3.0)
                s.connect((ESP_HOST, TCP_PORT_CMD))
                s.setsockopt(socket.IPPROTO_TCP, socket.TCP_NODELAY, 1)
                s.settimeout(2.0)
                self._cmd_sock = s
                return True
            except Exception as e:
                self.gain_status_label.setText(f"CMD conn err: {e}")
                return False

    def _send_cmd(self, cmd: str) -> str:
        """Send command and wait for response line. Returns response or empty string."""
        with self._cmd_lock:
            if self._cmd_sock is None:
                return ""
            try:
                self._cmd_sock.sendall((cmd + "\n").encode())
                # Read response (up to 256 bytes, expect \n terminated)
                data = b""
                while True:
                    chunk = self._cmd_sock.recv(256)
                    if not chunk:
                        break
                    data += chunk
                    if b"\n" in data:
                        break
                return data.decode().strip()
            except Exception as e:
                # Connection lost – reset
                try:
                    self._cmd_sock.close()
                except:
                    pass
                self._cmd_sock = None
                return f"ERR:{e}"

    def _on_gain_changed(self, index):
        """Called when user selects a new gain from the combo box."""
        gain = self.gain_combo.itemData(index)
        if gain is None:
            return
        self.gain_status_label.setText(f"Nastavuji GAIN x{gain} …")
        self.gain_status_label.setStyleSheet(
            "color: #ff0; font-family: Consolas, monospace; font-size: 12px; padding-left: 8px;"
        )
        # Run in background thread to avoid blocking GUI
        threading.Thread(target=self._apply_gain, args=(gain,), daemon=True).start()

    def _apply_gain(self, gain: int):
        """Send GAIN command to ESP32 (runs in background thread)."""
        if not self._connect_cmd():
            QtCore.QMetaObject.invokeMethod(
                self.gain_status_label, "setText",
                QtCore.Qt.QueuedConnection,
                QtCore.Q_ARG(str, "CMD spojeni selhalo")
            )
            return

        resp = self._send_cmd(f"GAIN:{gain}")
        if resp.startswith("GAIN_OK:"):
            confirmed_gain = int(resp.split(":")[1])
            # First update the LSB so new samples use the correct conversion
            self.receiver.current_gain = confirmed_gain
            self.receiver.lsb_mv = compute_lsb_mv(confirmed_gain)
            # Clear the ring buffer – old samples were converted with old gain LSB,
            # but samples captured between command send and actual HW change
            # would be decoded with wrong LSB → amplitude spike.
            # Short sleep lets the ESP32 flush any in-flight old-gain samples.
            time.sleep(0.1)  # ~50 samples @ 500 SPS – let old data drain
            self.ring.clear()
            # Update UI in main thread
            QtCore.QMetaObject.invokeMethod(
                self.gain_status_label, "setText",
                QtCore.Qt.QueuedConnection,
                QtCore.Q_ARG(str, f"GAIN x{confirmed_gain} potvrzeno ✓")
            )
            QtCore.QMetaObject.invokeMethod(
                self.gain_status_label, "setStyleSheet",
                QtCore.Qt.QueuedConnection,
                QtCore.Q_ARG(str, "color: #0f0; font-family: Consolas, monospace; font-size: 12px; padding-left: 8px;")
            )
        else:
            QtCore.QMetaObject.invokeMethod(
                self.gain_status_label, "setText",
                QtCore.Qt.QueuedConnection,
                QtCore.Q_ARG(str, f"GAIN chyba: {resp}")
            )
            QtCore.QMetaObject.invokeMethod(
                self.gain_status_label, "setStyleSheet",
                QtCore.Qt.QueuedConnection,
                QtCore.Q_ARG(str, "color: #f00; font-family: Consolas, monospace; font-size: 12px; padding-left: 8px;")
            )

    # ------------------------------------------------------------------ #
    def _on_rec_toggled(self, checked: bool):
        """Toggle CSV recording on/off."""
        if checked:
            self.recorder.start()
            fname = os.path.basename(self.recorder.filename)
            self.rec_status_label.setText(f"Nahrávám → {fname}")
            self.rec_status_label.setStyleSheet(
                "color: #f44; font-family: Consolas, monospace; font-size: 12px; padding-left: 8px;"
            )
            self.rec_btn.setText("⏹  STOP")
        else:
            n = self.recorder.sample_count
            fname = os.path.basename(self.recorder.filename)
            self.recorder.stop()
            self.rec_status_label.setText(f"Uloženo {n:,} vzorků → {fname} ✓")
            self.rec_status_label.setStyleSheet(
                "color: #0f0; font-family: Consolas, monospace; font-size: 12px; padding-left: 8px;"
            )
            self.rec_btn.setText("⏺  REC")

    def _update_rec_status(self):
        """Periodically update recording indicator."""
        if not self.recorder.recording:
            return
        n = self.recorder.sample_count
        elapsed = self.recorder.elapsed
        mins, secs = divmod(int(elapsed), 60)
        fname = os.path.basename(self.recorder.filename)
        self.rec_status_label.setText(
            f"⏺ {mins:02d}:{secs:02d}  |  {n:,} vzorků  |  {fname}"
        )

    # ------------------------------------------------------------------ #
    def _on_scene_click(self, event):
        """Handle double-click on a subplot to maximize/restore it."""
        if not event.double():
            return
        pos = event.scenePos()
        for i, p in enumerate(self.plots):
            if p.sceneBoundingRect().contains(pos):
                self._toggle_maximize(i)
                break

    def _toggle_maximize(self, index):
        """Toggle between maximized single-channel view and all-channels view."""
        if self._maximized_index == index:
            # --- Restore all plots ---
            self._maximized_index = None
            for i, p in enumerate(self.plots):
                p.setVisible(True)
                p.setMinimumHeight(40)
                p.setMaximumHeight(120)
                # Restore x-axis label visibility
                if i < NUM_CHANNELS - 1:
                    p.getAxis("bottom").setStyle(showValues=False)
                    p.getAxis("bottom").setLabel(None)
                else:
                    p.setLabel("bottom", "Čas", units="s")
                    p.getAxis("bottom").setStyle(showValues=True)
            # Re-link X axes
            for i in range(1, NUM_CHANNELS):
                self.plots[i].setXLink(self.plots[0])
        else:
            # --- Maximize selected plot ---
            self._maximized_index = index
            for i, p in enumerate(self.plots):
                if i == index:
                    p.setVisible(True)
                    p.setMinimumHeight(200)
                    p.setMaximumHeight(16777215)  # effectively unlimited
                    # Show x-axis on the maximized plot
                    p.setLabel("bottom", "Čas", units="s")
                    p.getAxis("bottom").setStyle(showValues=True)
                    # Thicker pen for maximized view
                    pen = pg.mkPen(color=CHANNEL_COLORS[i], width=2.0)
                    self.curves[i].setPen(pen)
                else:
                    p.setVisible(False)
                    p.setMinimumHeight(0)
                    p.setMaximumHeight(0)
                    # Restore normal pen width for hidden plots
                    pen = pg.mkPen(color=CHANNEL_COLORS[i], width=1.2)
                    self.curves[i].setPen(pen)

    # ------------------------------------------------------------------ #
    @staticmethod
    def _estimate_hr(signal: np.ndarray, fs: float) -> float:
        """Estimate heart rate (BPM) from a 1-D ECG signal using simple
        threshold-based R-peak detection.  Returns 0.0 if not enough data."""
        if signal.size < int(fs * 1.5):
            return 0.0

        # Band-pass-like preprocessing: remove DC, take abs of derivative
        sig = signal.astype(np.float64)
        sig -= np.mean(sig)
        diff = np.diff(sig)
        energy = diff * diff  # squared first-difference → accentuates QRS

        # Adaptive threshold: 40 % of the 95th-percentile energy
        thr = np.percentile(energy, 95) * 0.40
        if thr <= 0:
            return 0.0

        # Min distance between R-peaks: 300 ms (= 200 BPM max)
        min_dist = int(fs * 0.30)

        # Find peaks above threshold with minimum distance
        above = energy > thr
        peaks = []
        i = 0
        n = len(above)
        while i < n:
            if above[i]:
                # Find local maximum in this supra-threshold region
                best = i
                while i < n and above[i]:
                    if energy[i] > energy[best]:
                        best = i
                    i += 1
                if not peaks or (best - peaks[-1]) >= min_dist:
                    peaks.append(best)
            else:
                i += 1

        if len(peaks) < 2:
            return 0.0

        # Average RR interval → BPM
        rr_samples = np.diff(peaks).astype(np.float64)
        rr_sec = np.mean(rr_samples) / fs
        bpm = 60.0 / rr_sec
        # Sanity clamp
        if bpm < 20 or bpm > 250:
            return 0.0
        return bpm

    # ------------------------------------------------------------------ #
    def _on_pause_toggled(self):
        """Zastavit nebo znovu spustit vykreslování grafu."""
        self._plot_paused = self.pause_btn.isChecked()
        if self._plot_paused:
            self.timer.stop()
            self.pause_btn.setText("▶ Spustit vykreslování")
        else:
            self.timer.start(PLOT_INTERVAL_MS)
            self.pause_btn.setText("⏸ Pozastavit vykreslování")

    def _update(self):
        """Called every PLOT_INTERVAL_MS – update curves & status."""
        if self._plot_paused:
            return
        t, d, ts_us = self.ring.snapshot()
        if t.size == 0:
            self._x_view_min = self._x_view_max = None  # reset pro plynulý start po opětovném příjmu
            self._update_status(0, np.array([], dtype=np.uint64), None)
            return

        # Sanitize: nahradit nan/inf konečnými hodnotami, aby PyQtGraph ViewBox nezpůsobil "overflow in cast"
        t = np.nan_to_num(t.astype(np.float64), nan=0.0, posinf=0.0, neginf=0.0).astype(np.float32)
        for i in range(NUM_CHANNELS):
            ch = np.asarray(d[i], dtype=np.float64)
            ch = np.nan_to_num(ch, nan=0.0, posinf=0.0, neginf=0.0).astype(np.float32)
            self.curves[i].setData(t, ch)
            # # --- Y_SPAN: zakomentováno – aby fungovalo manuální škálování kolečkem myši ---
            # # Rozsah Y: buď z Y_SPAN_PER_CHANNEL_MV[i], nebo auto z dat (min. Y_SPAN_MIN_MV)
            # if ch.size > 0:
            #     lo, hi = float(np.min(ch)), float(np.max(ch))
            #     center = (lo + hi) * 0.5
            #     data_span = hi - lo
            #     if Y_SPAN_PER_CHANNEL_MV is not None and len(Y_SPAN_PER_CHANNEL_MV) > i:
            #         span = float(Y_SPAN_PER_CHANNEL_MV[i])
            #     else:
            #         span = max(data_span, Y_SPAN_MIN_MV)
            #     half = span * 0.5
            #     pad = span * 0.05
            #     self.plots[i].setYRange(center - half - pad, center + half + pad, padding=0)

        # Plynulý scroll osy X: plynulé přibližování k cílovému rozsahu
        t_newest = float(t[-1])
        t_oldest = t_newest - WINDOW_SEC
        if not np.isfinite(t_newest):
            t_newest = 0.0
            t_oldest = -WINDOW_SEC
        if not np.isfinite(t_oldest):
            t_oldest = t_newest - WINDOW_SEC
        alpha = SCROLL_SMOOTH_ALPHA
        if self._x_view_min is None or self._x_view_max is None:
            self._x_view_min, self._x_view_max = t_oldest, t_newest
        else:
            self._x_view_min += alpha * (t_oldest - self._x_view_min)
            self._x_view_max += alpha * (t_newest - self._x_view_max)
        for p in self.plots:
            p.setXRange(self._x_view_min, self._x_view_max, padding=0)

        # Pass full-res Lead I (channel 0) for HR estimation
        self._update_status(t.size, ts_us, d[0] if d.shape[1] > 0 else None)

    def _update_status(self, n_samples, ts_us, lead_i_data):
        now = time.time()
        dt = now - self._last_stat_time
        pkts = self.receiver.packets_received
        if dt >= 0.5:
            self._last_rate = (pkts - self._last_pkt_count) / dt
            self._last_stat_time = now
            self._last_pkt_count = pkts
        rate = getattr(self, '_last_rate', 0.0)

        # Compute inter-sample jitter from ESP32 timestamps
        jitter_str = "—"
        if ts_us.size > 10:
            diffs = np.diff(ts_us.astype(np.int64))
            if diffs.size > 0:
                mean_dt_us = np.mean(diffs)
                std_dt_us  = np.std(diffs)
                jitter_str = f"Δt={mean_dt_us:.0f}±{std_dt_us:.0f} µs"

        latency_ms = self.receiver.last_latency_us / 1000.0

        # Heart rate from Lead I (throttled to ~1 Hz – expensive)
        if now - self._last_hr_time >= 1.0:
            self._last_hr_time = now
            if lead_i_data is not None and lead_i_data.size > 0:
                bpm = self._estimate_hr(lead_i_data, SAMPLE_RATE)
                if bpm > 0:
                    self._last_hr_str = f"HR: {bpm:.0f} BPM"
                else:
                    self._last_hr_str = "HR: —"
        hr_str = self._last_hr_str

        conn_icon = "🟢" if self.receiver.connected else "🔴"
        gain_str = f"Gain:x{self.receiver.current_gain}"
        self.status_label.setText(
            f"{conn_icon}  {self.receiver.status}   |   "
            f"Vzorků: {pkts:,}   |   "
            f"Frames: {self.receiver.frames_received}   |   "
            f"{rate:.0f} smp/s   |   "
            f"Buf: {n_samples}/{WINDOW_SAMPLES}   |   "
            f"Lat: {latency_ms:.1f}ms   |   "
            f"{jitter_str}   |   "
            f"{gain_str}   |   "
            f"{hr_str}"
        )

    def closeEvent(self, event):
        self.receiver.stop()
        if self.recorder.recording:
            self.recorder.stop()
        event.accept()


# =============================================================================
# MAIN
# =============================================================================

def main():
    # Create ring buffer
    ring = RingBuffer(WINDOW_SAMPLES, NUM_CHANNELS)

    # Create CSV recorder (shared between receiver & GUI)
    recorder = CsvRecorder()

    # Start TCP receiver thread
    receiver = TcpReceiver(ring, recorder)
    receiver.start()

    # Launch Qt application
    app = QtWidgets.QApplication(sys.argv)
    app.setStyle("Fusion")

    # Dark palette
    palette = QtGui.QPalette()
    palette.setColor(QtGui.QPalette.Window,          QtGui.QColor(30, 30, 30))
    palette.setColor(QtGui.QPalette.WindowText,      QtGui.QColor(200, 200, 200))
    palette.setColor(QtGui.QPalette.Base,            QtGui.QColor(20, 20, 20))
    palette.setColor(QtGui.QPalette.AlternateBase,   QtGui.QColor(40, 40, 40))
    palette.setColor(QtGui.QPalette.Text,            QtGui.QColor(200, 200, 200))
    palette.setColor(QtGui.QPalette.Button,          QtGui.QColor(50, 50, 50))
    palette.setColor(QtGui.QPalette.ButtonText,      QtGui.QColor(200, 200, 200))
    app.setPalette(palette)

    win = ECGPlotWindow(ring, receiver, recorder)
    win.show()

    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
