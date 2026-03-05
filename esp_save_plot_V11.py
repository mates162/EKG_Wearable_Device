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
    pip install pyqtgraph PyQt5 numpy scipy
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

try:
    from scipy.signal import butter, filtfilt, iirnotch, sosfiltfilt
    _SCIPY_AVAILABLE = True
except ImportError:
    _SCIPY_AVAILABLE = False

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

# --- Y osa: fixní rozptyl -1..1 mV jen při zapnutých filtrech; při vypnutých filtrech autoscale ---
Y_SPAN_MV = 2.0   # rozptyl v mV (rozsah -1 .. 1 mV)
Y_VIEW_MIN_MV = -1.0
Y_VIEW_MAX_MV = 1.0
# Jak často obnovit fixní rozsah Y (s), aby to neškubalo a nezpomalovalo
Y_FIXED_RANGE_INTERVAL_S = 0.4

# TCP receive buffer size
TCP_RECV_SIZE   = 16384

# ---- Binary protocol constants (must match ESP32 firmware v11) ----
SYNC_MARKER     = bytes([0xA5, 0x5A])
RAW_SPI_BYTES   = 27          # 3B status + 8×3B channels
SAMPLE_SIZE     = 4 + RAW_SPI_BYTES   # 31 bytes: 4B ts + 27B raw
FRAME_HEADER    = 3           # 2B sync + 1B count
SAMPLES_PER_FRAME = 10
FRAME_EXPECTED  = FRAME_HEADER + SAMPLES_PER_FRAME * SAMPLE_SIZE  # 623

# ---- ADS1298 → Voltage conversion (vše v mV) ----
# V_input = ADC_code × VREF / (Gain × 2^23); výstup převádíme na mV (× 1e3)
ADS_VREF        = 2.4              # Internal reference voltage [V]
ADS_GAIN        = 6                # PGA gain (CHnSET[6:4] = 000)
ADS_LSB_MV      = (ADS_VREF / (ADS_GAIN * (2**23))) * 1e3   # mV per LSB

# Valid ADS1298 PGA gain values
VALID_GAINS     = [1, 2, 3, 4, 6, 8, 12]

# ---- Klinické EKG filtry (AAMI / diagnostické EKG) ----
# Pásmový filtr: odstraní baseline wander (<0.5 Hz) a vysokofrekvenční šum (>40 Hz)
FILTER_BANDPASS_LOW   = 0.5   # Hz
FILTER_BANDPASS_HIGH  = 40.0  # Hz
FILTER_BANDPASS_ORDER = 4     # řád Butterworth
# Notch: potlačení síťového rušení (50 Hz EU / 60 Hz US)
NOTCH_QUALITY         = 30    # Q pro úzký zářez
# Minimální počet vzorků pro stabilní filtrování (filtfilt)
FILTER_MIN_SAMPLES    = int(SAMPLE_RATE * 1.0)  # alespoň 1 s

def compute_lsb_mv(gain: int) -> float:
    """Vrátí mV na LSB pro dané PGA gain (vstup i výstup zůstávají v mV)."""
    return (ADS_VREF / (gain * (2**23))) * 1e3


# =============================================================================
# KLINICKÉ EKG FILTRY
# =============================================================================

def apply_ecg_filters(
    data: np.ndarray,
    fs: float,
    bandpass_low: float = FILTER_BANDPASS_LOW,
    bandpass_high: float = FILTER_BANDPASS_HIGH,
    notch_hz: float | None = 50.0,
) -> np.ndarray:
    """Aplikuje filtry používané při záznamu klinického EKG:
    - Pásmový filtr 0,5–40 Hz (Butterworth): baseline wander + HF šum
    - Notch 50/60 Hz: síťové rušení

    data: shape (num_channels, N) v mV
    fs: vzorkovací kmitočet [Hz]
    notch_hz: None = bez notch, jinak 50 nebo 60
    Vrací filtrovaná data stejného tvaru. Při nedostupném scipy vrací data beze změny.
    """
    if not _SCIPY_AVAILABLE:
        return data
    if data.size == 0:
        return data
    num_ch, n = data.shape
    if n < FILTER_MIN_SAMPLES:
        return data

    out = np.empty_like(data, dtype=np.float64)
    for ch in range(num_ch):
        sig = np.asarray(data[ch], dtype=np.float64)

        # 1) Bandpass 0.5–40 Hz (AAMI EC11)
        nyq = 0.5 * fs
        low = max(0.01, bandpass_low / nyq)
        high = min(0.99, bandpass_high / nyq)
        if low < high:
            sos = butter(FILTER_BANDPASS_ORDER, [low, high], btype="band", output="sos")
            sig = sosfiltfilt(sos, sig, axis=-1)

        # 2) Notch 50/60 Hz
        if notch_hz is not None and 0 < notch_hz < nyq:
            w0 = notch_hz / nyq
            if w0 < 0.99 and w0 > 0.01:
                b, a = iirnotch(notch_hz, NOTCH_QUALITY, fs)
                sig = filtfilt(b, a, sig, axis=-1)

        out[ch] = sig

    return out.astype(data.dtype, copy=False)


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

# Sloupce s kanály (pro načtení CSV) – pořadí musí odpovídat 12 svodům
CSV_CHANNEL_COLUMNS = [
    "I_mV", "II_mV", "III_mV", "aVR_mV", "aVL_mV", "aVF_mV",
    "V1_mV", "V2_mV", "V3_mV", "V4_mV", "V5_mV", "V6_mV",
]


def load_csv(filepath: str):
    """Načte CSV soubor z aplikace (timestamp_us + 12 kanálů mV).
    Vrací (time_sec, data) kde time_sec je 1D pole v sekundách od začátku,
    data má tvar (12, N) v mV. Při chybě vrací None."""
    try:
        with open(filepath, "r", newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            fieldnames = list(reader.fieldnames or [])
            rows = list(reader)
        if not rows or not fieldnames:
            return None

        # Sloupec s časem
        ts_key = "timestamp_us" if "timestamp_us" in fieldnames else None
        if ts_key is None:
            for key in fieldnames:
                if "timestamp" in key.lower():
                    ts_key = key
                    break
        if ts_key is None:
            return None

        # Sloupce kanálů (12 svodů v pořadí I, II, III, aVR, aVL, aVF, V1–V6)
        ch_cols = [c for c in CSV_CHANNEL_COLUMNS if c in fieldnames]
        if len(ch_cols) != NUM_CHANNELS:
            ch_cols = [c for c in fieldnames if c and c not in ("timestamp_us", "gain")][:NUM_CHANNELS]
        if len(ch_cols) != NUM_CHANNELS:
            return None

        ts_list = []
        ch_arrays = [[] for _ in range(NUM_CHANNELS)]
        for row in rows:
            try:
                t = int(float(row.get(ts_key, 0)))
            except (ValueError, TypeError):
                continue
            ts_list.append(t)
            for i, col in enumerate(ch_cols):
                try:
                    ch_arrays[i].append(float(row.get(col, 0)))
                except (ValueError, TypeError):
                    ch_arrays[i].append(0.0)

        if not ts_list:
            return None
        ts_us = np.array(ts_list, dtype=np.uint64)
        time_sec = (ts_us.astype(np.float64) - float(ts_us[0])) / 1e6
        data = np.array(ch_arrays, dtype=np.float32)  # (12, N)
        return time_sec, data
    except Exception:
        return None


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
            # Jednotky "mV" v textu popisku – bez units= v API, aby PyQtGraph neaplikoval SI předpony (mmV) a neškáloval hodnoty
            p.setLabel("left", f"{CHANNEL_LABELS[i]}, mV", units="")
            left_axis = p.getAxis("left")
            left_axis.setStyle(tickTextWidth=52, autoExpandTextSpace=False)
            # Zakázat automatické SI předpony na ose Y – hodnoty zůstávají v mV (0.2, 0.5, 1.0), ne ve stovkách
            if hasattr(left_axis, "enableAutoSIPrefix"):
                left_axis.enableAutoSIPrefix(False)
            p.showGrid(x=True, y=True, alpha=0.25)
            p.setMouseEnabled(x=False, y=True)
            p.enableAutoRange(axis="y", enable=False)
            p.setYRange(Y_VIEW_MIN_MV, Y_VIEW_MAX_MV, padding=0)
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

        # --- Y: při zapnutých filtrech fixní -1..1 mV (kolečko pustí, Autoscale vrátí); při vypnutých filtrech autoscale ---
        self._y_autoscale = True
        self._applying_autoscale = False
        self._last_y_fixed_range_time = 0.0
        for p in self.plots:
            p.getViewBox().sigRangeChanged.connect(self._on_y_range_changed)

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

        # --- Klinické filtry EKG (0,5–40 Hz + notch) ---
        controls_layout.addSpacing(24)
        self.filter_cb = QtWidgets.QCheckBox("Klinické filtry: ZAPNUTO")
        self.filter_cb.setChecked(True)
        self.filter_cb.setStyleSheet(
            "QCheckBox { font-size: 13px; font-weight: bold; }"
            "QCheckBox { color: #0c5; }"
            "QCheckBox:unchecked { color: #666; }"
            "QCheckBox::indicator { width: 20px; height: 20px; border: 2px solid #555; border-radius: 3px; background: #2a2a2a; }"
            "QCheckBox:checked::indicator { background: #0a5; border-color: #0f8; }"
            "QCheckBox:unchecked::indicator { background: #333; }"
        )
        self.filter_cb.setToolTip(
            "Pásmový filtr 0,5–40 Hz (baseline + HF šum) a notch 50/60 Hz (síť). Doporučeno pro diagnostické EKG."
        )
        self.filter_cb.stateChanged.connect(self._on_filter_toggled)
        controls_layout.addWidget(self.filter_cb)
        self.notch_combo = QtWidgets.QComboBox()
        self.notch_combo.setFixedWidth(70)
        self.notch_combo.addItem("50 Hz", 50.0)
        self.notch_combo.addItem("60 Hz", 60.0)
        self.notch_combo.setCurrentIndex(0)
        self.notch_combo.setStyleSheet(
            "QComboBox { background: #333; color: #0af; font-size: 12px; border: 1px solid #555; padding: 2px 4px; }"
        )
        self.notch_combo.setToolTip("Frekvence notch filtru (síťové rušení)")
        controls_layout.addWidget(QtWidgets.QLabel("Notch:"))
        controls_layout.addWidget(self.notch_combo)

        # --- Autoscale Y: návrat na rozptyl -1 .. 1 mV ---
        controls_layout.addSpacing(24)
        self.autoscale_btn = QtWidgets.QPushButton("Autoscale Y")
        self.autoscale_btn.setFixedWidth(110)
        self.autoscale_btn.setStyleSheet(
            "QPushButton { background: #333; color: #0af; font-size: 13px; "
            "font-weight: bold; border: 1px solid #555; padding: 2px 8px; }"
        )
        self.autoscale_btn.setToolTip("Nastaví osu Y všech kanálů na -1 až 1 mV (rozptyl 2 mV)")
        self.autoscale_btn.clicked.connect(self._on_autoscale_clicked)
        controls_layout.addWidget(self.autoscale_btn)

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

        # --- Načtený CSV: zobrazení a procházení ---
        self._csv_data = None   # (time_sec: np.ndarray, data: (12, N)) nebo None
        self._csv_path = None
        self._csv_window_start = 0.0   # začátek okna v sekundách (vztaženo na začátek souboru)
        self._view_mode = "live"       # "live" | "csv"

        # --- Tlačítko Načíst CSV ---
        self.load_csv_btn = QtWidgets.QPushButton("📂 Načíst CSV")
        self.load_csv_btn.setFixedWidth(120)
        self.load_csv_btn.setStyleSheet(
            "QPushButton { background: #333; color: #ccc; font-size: 13px; "
            "font-weight: bold; border: 1px solid #555; padding: 2px 8px; }"
        )
        self.load_csv_btn.clicked.connect(self._on_load_csv)
        controls_layout.addWidget(self.load_csv_btn)

        # --- Režim zobrazení: Živý / CSV ---
        self.view_mode_combo = QtWidgets.QComboBox()
        self.view_mode_combo.setFixedWidth(140)
        self.view_mode_combo.addItem("Živý signál", "live")
        self.view_mode_combo.addItem("Načtený CSV", "csv")
        self.view_mode_combo.setCurrentIndex(0)
        self.view_mode_combo.setEnabled(False)
        self.view_mode_combo.setStyleSheet(
            "QComboBox { background: #333; color: #0af; font-size: 12px; "
            "border: 1px solid #555; padding: 2px 6px; }"
        )
        self.view_mode_combo.currentIndexChanged.connect(self._on_view_mode_changed)
        controls_layout.addWidget(QtWidgets.QLabel("Režim:"))
        controls_layout.addWidget(self.view_mode_combo)

        # --- Slider procházení CSV (viditelný jen v režimu CSV) ---
        self.csv_slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.csv_slider.setMinimum(0)
        self.csv_slider.setMaximum(0)
        self.csv_slider.setValue(0)
        self.csv_slider.setMinimumWidth(200)
        self.csv_slider.valueChanged.connect(self._on_csv_slider_changed)
        self.csv_slider.setVisible(False)
        controls_layout.addWidget(self.csv_slider, stretch=0)

        self.csv_pos_label = QtWidgets.QLabel("")
        self.csv_pos_label.setStyleSheet("color: #888; font-size: 11px;")
        self.csv_pos_label.setVisible(False)
        controls_layout.addWidget(self.csv_pos_label)

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
    # Channels to try for HR auto-selection ((0)I, (1)II, (2)III, (5)aVF – limb leads)
    HR_CANDIDATE_CHANNELS = (0, 1)

    @classmethod
    def _estimate_hr(cls, signal: np.ndarray, fs: float) -> float:
        """Estimate heart rate (BPM) from a 1-D ECG signal using simple
        threshold-based R-peak detection.  Returns 0.0 if not enough data."""
        bpm, _ = cls._estimate_hr_with_quality(signal, fs)
        return bpm

    @classmethod
    def _estimate_hr_with_quality(cls, signal: np.ndarray, fs: float) -> tuple:
        """Estimate HR and quality (stability of RR intervals).
        Returns (bpm, quality). quality is in [0, 1], higher = more stable RR.
        Returns (0.0, 0.0) if not enough data or invalid."""
        if signal.size < int(fs * 1.5):
            return (0.0, 0.0)

        sig = signal.astype(np.float64)
        sig -= np.mean(sig)
        diff = np.diff(sig)
        energy = diff * diff

        thr = np.percentile(energy, 95) * 0.40
        if thr <= 0:
            return (0.0, 0.0)

        min_dist = int(fs * 0.30)
        above = energy > thr
        peaks = []
        i = 0
        n = len(above)
        while i < n:
            if above[i]:
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
            return (0.0, 0.0)

        rr_samples = np.diff(peaks).astype(np.float64)
        rr_sec = np.mean(rr_samples) / fs
        bpm = 60.0 / rr_sec
        if bpm < 20 or bpm > 250:
            return (0.0, 0.0)
        # Quality: inverse of coefficient of variation of RR (stable = high quality)
        rr_std = np.std(rr_samples)
        rr_mean = np.mean(rr_samples)
        if rr_mean <= 0:
            return (bpm, 0.0)
        cv = rr_std / rr_mean
        quality = 1.0 / (1.0 + cv)
        return (bpm, quality)

    @classmethod
    def _best_hr_from_channels(cls, data: np.ndarray, fs: float) -> tuple:
        """Try limb leads (I, II, III, aVF) and return (bpm, channel_index) with best quality.
        data shape: (num_channels, N). Returns (0.0, 0) if no valid HR."""
        if data is None or data.ndim != 2 or data.shape[1] < int(fs * 1.5):
            return (0.0, 0)
        best_bpm = 0.0
        best_quality = 0.0
        best_ch = 0
        for ch in cls.HR_CANDIDATE_CHANNELS:
            if ch >= data.shape[0]:
                continue
            bpm, quality = cls._estimate_hr_with_quality(data[ch], fs)
            if bpm > 0 and quality > best_quality:
                best_quality = quality
                best_bpm = bpm
                best_ch = ch
        return (best_bpm, best_ch)

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

    def _on_filter_toggled(self, _state):
        """Při zapnutí filtrů: fixní Y -1..1 mV; při vypnutí: Y autoscale. Aktualizovat text checkboxu."""
        filters_on = self.filter_cb.isChecked()
        if filters_on:
            self.filter_cb.setText("Klinické filtry: ZAPNUTO")
            self.autoscale_btn.setEnabled(True)
            self.autoscale_btn.setToolTip("Nastaví osu Y všech kanálů na -1 až 1 mV (rozptyl 2 mV)")
            for p in self.plots:
                p.enableAutoRange(axis="y", enable=False)
            self._y_autoscale = True
            self._applying_autoscale = True
            try:
                for p in self.plots:
                    p.setYRange(Y_VIEW_MIN_MV, Y_VIEW_MAX_MV, padding=0)
            finally:
                self._applying_autoscale = False
            self._last_y_fixed_range_time = time.time()
        else:
            self.filter_cb.setText("Klinické filtry: VYPNUTO")
            self.autoscale_btn.setEnabled(False)
            self.autoscale_btn.setToolTip("Dostupné jen při zapnutých klinických filtrech (Y má pak autoscale)")
            for p in self.plots:
                p.enableAutoRange(axis="y", enable=True)
            self._y_autoscale = False

    def _on_y_range_changed(self, vb, _range):
        """Když uživatel změní rozsah Y (kolečko/drag) a filtry jsou zapnuté, vypnout fixní -1..1 mV."""
        if not self.filter_cb.isChecked() or self._applying_autoscale:
            return
        yr = vb.viewRange()[1]
        ymin, ymax = yr[0], yr[1]
        tol = 0.05
        if abs(ymin - Y_VIEW_MIN_MV) > tol or abs(ymax - Y_VIEW_MAX_MV) > tol:
            self._y_autoscale = False

    def _on_autoscale_clicked(self):
        """Nastavit všechny kanály na rozptyl -1 .. 1 mV."""
        self._y_autoscale = True
        self._applying_autoscale = True
        try:
            for p in self.plots:
                p.setYRange(Y_VIEW_MIN_MV, Y_VIEW_MAX_MV, padding=0)
        finally:
            self._applying_autoscale = False

    def _on_load_csv(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "Načíst CSV", LOGS_DIR or SCRIPT_DIR,
            "CSV soubory (*.csv);;Všechny soubory (*)"
        )
        if not path:
            return
        result = load_csv(path)
        if result is None:
            QtWidgets.QMessageBox.warning(
                self, "Chyba", "Soubor se nepodařilo načíst nebo nemá očekávaný formát\n(timestamp_us + 12 kanálů mV)."
            )
            return
        time_sec, data = result
        self._csv_data = (time_sec, data)
        self._csv_path = path
        self._csv_window_start = 0.0
        self.setWindowTitle(f"ESP32 12-Lead ECG — {os.path.basename(path)}")
        self.view_mode_combo.setEnabled(True)
        self.view_mode_combo.setCurrentIndex(1)  # přepnout na "Načtený CSV"
        self._view_mode = "csv"
        total_sec = float(time_sec[-1]) if time_sec.size else 0.0
        max_start = max(0, total_sec - WINDOW_SEC)
        self.csv_slider.setMaximum(int(max_start * 10))  # krok 0.1 s
        self.csv_slider.setValue(0)
        self.csv_slider.setVisible(True)
        self.csv_pos_label.setVisible(True)
        self._update_csv_pos_label()

    def _on_view_mode_changed(self, index):
        mode = self.view_mode_combo.itemData(index)
        self._view_mode = mode if mode else "live"
        self.csv_slider.setVisible(self._view_mode == "csv" and self._csv_data is not None)
        self.csv_pos_label.setVisible(self._view_mode == "csv" and self._csv_data is not None)
        if self._view_mode == "csv":
            self._update_csv_pos_label()
            if self._csv_path:
                self.setWindowTitle(f"ESP32 12-Lead ECG — {os.path.basename(self._csv_path)}")
        else:
            self.setWindowTitle("ESP32 12-Lead ECG  —  Real-Time WiFi Plotter (v11)")

    def _on_csv_slider_changed(self, value):
        self._csv_window_start = value / 10.0
        self._update_csv_pos_label()

    def _update_csv_pos_label(self):
        if self._csv_data is None:
            self.csv_pos_label.setText("")
            return
        time_sec, data = self._csv_data
        n = time_sec.size
        total_sec = float(time_sec[-1]) if n else 0.0
        self.csv_pos_label.setText(
            f"Čas: {self._csv_window_start:.1f}–{min(self._csv_window_start + WINDOW_SEC, total_sec):.1f} s / {total_sec:.1f} s  ({n} vzorků)"
        )

    def _update(self):
        """Called every PLOT_INTERVAL_MS – update curves & status."""
        if self._plot_paused:
            return

        # Režim načteného CSV: zobrazit výřez podle slideru
        if self._view_mode == "csv" and self._csv_data is not None:
            time_sec, data = self._csv_data
            t_end = self._csv_window_start + WINDOW_SEC
            mask = (time_sec >= self._csv_window_start) & (time_sec <= t_end)
            t = time_sec[mask].astype(np.float32)
            if t.size == 0:
                self._update_status(0, np.array([], dtype=np.uint64), None)
                return
            d = data[:, mask].astype(np.float64)
            # Klinické filtry EKG při zobrazení CSV
            if self.filter_cb.isChecked():
                notch_hz = self.notch_combo.currentData()
                d = apply_ecg_filters(d, SAMPLE_RATE, notch_hz=notch_hz)
            t = np.nan_to_num(t.astype(np.float64), nan=0.0, posinf=0.0, neginf=0.0).astype(np.float32)
            for i in range(NUM_CHANNELS):
                ch = np.asarray(d[i], dtype=np.float64)
                ch = np.nan_to_num(ch, nan=0.0, posinf=0.0, neginf=0.0).astype(np.float32)
                self.curves[i].setData(t, ch)
            if self.filter_cb.isChecked() and self._y_autoscale:
                now = time.time()
                if now - self._last_y_fixed_range_time >= Y_FIXED_RANGE_INTERVAL_S:
                    self._last_y_fixed_range_time = now
                    self._applying_autoscale = True
                    try:
                        for p in self.plots:
                            p.setYRange(Y_VIEW_MIN_MV, Y_VIEW_MAX_MV, padding=0)
                    finally:
                        self._applying_autoscale = False
            for p in self.plots:
                p.setXRange(self._csv_window_start, min(t_end, float(time_sec[-1])), padding=0)
            self._update_status(d.shape[1], np.array([], dtype=np.uint64), d if d.shape[1] > 0 else None)
            return

        t, d, ts_us = self.ring.snapshot()
        if t.size == 0:
            self._x_view_min = self._x_view_max = None  # reset pro plynulý start po opětovném příjmu
            self._update_status(0, np.array([], dtype=np.uint64), None)
            return

        # Klinické filtry EKG (0,5–40 Hz + notch)
        if self.filter_cb.isChecked():
            notch_hz = self.notch_combo.currentData()
            d = apply_ecg_filters(d, SAMPLE_RATE, notch_hz=notch_hz)

        # Sanitize: nahradit nan/inf konečnými hodnotami, aby PyQtGraph ViewBox nezpůsobil "overflow in cast"
        t = np.nan_to_num(t.astype(np.float64), nan=0.0, posinf=0.0, neginf=0.0).astype(np.float32)
        for i in range(NUM_CHANNELS):
            ch = np.asarray(d[i], dtype=np.float64)
            ch = np.nan_to_num(ch, nan=0.0, posinf=0.0, neginf=0.0).astype(np.float32)
            self.curves[i].setData(t, ch)

        # Y: při zapnutých filtrech a fixním režimu držet -1..1 mV, obnovovat jen občas (ne každý snímek)
        if self.filter_cb.isChecked() and self._y_autoscale:
            now = time.time()
            if now - self._last_y_fixed_range_time >= Y_FIXED_RANGE_INTERVAL_S:
                self._last_y_fixed_range_time = now
                self._applying_autoscale = True
                try:
                    for p in self.plots:
                        p.setYRange(Y_VIEW_MIN_MV, Y_VIEW_MAX_MV, padding=0)
                finally:
                    self._applying_autoscale = False

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

        # Pass full-res 12-channel data for HR auto-selection (best of I, II, III, aVF)
        self._update_status(t.size, ts_us, d if d.shape[1] > 0 else None)

    def _update_status(self, n_samples, ts_us, ecg_data):
        """ecg_data: None, 1D array (single channel), or 2D (num_channels, N) for HR auto-selection."""
        if self._view_mode == "csv" and self._csv_data is not None:
            time_sec, data = self._csv_data
            total_sec = float(time_sec[-1]) if time_sec.size else 0.0
            hr_str = "HR: —"
            if ecg_data is not None and np.size(ecg_data) > 0:
                if ecg_data.ndim == 2 and ecg_data.shape[0] >= 6:
                    bpm, ch = self._best_hr_from_channels(ecg_data, SAMPLE_RATE)
                else:
                    sig = ecg_data.ravel()
                    bpm = self._estimate_hr(sig, SAMPLE_RATE)
                    ch = 0
                if bpm > 0:
                    lead_name = ("I", "II", "III", "aVR", "aVL", "aVF")[ch] if ch < 6 else str(ch)
                    hr_str = f"HR: {bpm:.0f} BPM ({lead_name})"
            self.status_label.setText(
                f"📂 CSV: {os.path.basename(self._csv_path or '')}   |   "
                f"Čas: {self._csv_window_start:.1f}–{min(self._csv_window_start + WINDOW_SEC, total_sec):.1f} s / {total_sec:.1f} s   |   "
                f"Vzorků: {n_samples:,}   |   {hr_str}"
            )
            return
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

        # Heart rate with auto-selected best channel (I, II, III, aVF); throttled ~1 Hz
        if now - self._last_hr_time >= 1.0:
            self._last_hr_time = now
            if ecg_data is not None and np.size(ecg_data) > 0:
                if ecg_data.ndim == 2 and ecg_data.shape[0] >= 6:
                    bpm, ch = self._best_hr_from_channels(ecg_data, SAMPLE_RATE)
                else:
                    sig = ecg_data.ravel()
                    bpm = self._estimate_hr(sig, SAMPLE_RATE)
                    ch = 0
                if bpm > 0:
                    lead_name = ("I", "II", "III", "aVR", "aVL", "aVF")[ch] if ch < 6 else str(ch)
                    self._last_hr_str = f"HR: {bpm:.0f} BPM ({lead_name})"
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
