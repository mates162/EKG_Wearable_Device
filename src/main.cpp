/*
 * ESP32-C6  12-Lead ECG  (ADS1298 + MPU9250 disabled)  v10 – spike fix
 * --------------------------------------------------------------------
 * Sampling rate : 500 Hz  (ADS1298 HR mode, DR=110)
 * Serial baud   : 250 000  (debug only)
 * WiFi AP mode  : SSID "ESP32_ECG", TCP port 3333 (data→PC), 3334 (PC→ESP, reserved)
 * SPI bus shared between ADS1298 (CS4) and MPU9250 (CS3, kept HIGH)
 *
 * ADS1298 electrode mapping (12-lead ECG):
 *   CH1: IN1P=V1   IN1N=WCT
 *   CH2: IN2P=LA   IN2N=RA     → Lead I
 *   CH3: IN3P=LL   IN3N=RA     → Lead II
 *   CH4: IN4P=V2   IN4N=WCT
 *   CH5: IN5P=V3   IN5N=WCT
 *   CH6: IN6P=V4   IN6N=WCT
 *   CH7: IN7P=V5   IN7N=WCT
 *   CH8: IN8P=V6   IN8N=WCT
 *
 * Internal WCT = (RA + LA + LL) / 3 via WCTA/B/C amplifiers.
 * Remaining leads derived in software:
 *   III = II − I,  aVR = −(I+II)/2,  aVL = I − II/2,  aVF = II − I/2
 *
 * v10 changes (spike elimination):
 *   - Atomic SPI block transfer (SPI.transferBytes) instead of byte-by-byte
 *   - Interrupts disabled during SPI read to prevent WiFi IRQ mid-transfer
 *   - DRDY ISR captures timestamp at exact falling edge
 *   - Dedicated high-priority FreeRTOS task for SPI reads
 *   - Separate lower-priority task for TCP frame assembly & send
 *   - Ring buffer (64 samples) decouples SPI read from TCP send
 *   - ADS1298 status byte validation (0xCx) – bad samples discarded
 *   - Overflow detection & debug stats
 *
 * v11 changes (RDATA mode):
 *   - Replaced RDATAC (continuous read) with explicit RDATA command per sample
 *   - Each read: DRDY fires → RDATA opcode sent → 27 bytes clocked out
 *   - More stable: no risk of read desynchronisation in continuous mode
 *   - Simplified gain change: no SDATAC/RDATAC transitions needed
 *   - Register access possible at any time (device not in RDATAC lock)
 */

#include <Arduino.h>
#include <SPI.h>
#include <WiFi.h>
#include <WiFiAP.h>
#include <esp_sleep.h>
#include <freertos/FreeRTOS.h>
#include <freertos/task.h>
#include <freertos/semphr.h>

// ======================== PIN DEFINITIONS ========================
#define PIN_SPI_MISO     2
#define PIN_SPI_MOSI     7
#define PIN_SPI_CLK      6
#define PIN_CS_ADS      20   // SPI_CS4 → ADS1298
#define PIN_CS_MPU      19   // SPI_CS3 → MPU9250  (kept HIGH = disabled)
#define PIN_ADS_DRDY    21
#define PIN_ADS_PWDN    18
#define PIN_ADS_RESET   22
#define PIN_ADS_START   23
#define PIN_MPU_INT      3

// ======================== ADS1298 SPI OPCODES ====================
#define ADS_CMD_WAKEUP   0x02
#define ADS_CMD_STANDBY  0x04
#define ADS_CMD_RESET    0x06
#define ADS_CMD_START    0x08
#define ADS_CMD_STOP     0x0A
#define ADS_CMD_RDATAC   0x10
#define ADS_CMD_SDATAC   0x11
#define ADS_CMD_RDATA    0x12

// ======================== ADS1298 REGISTER MAP ===================
#define ADS_REG_ID          0x00
#define ADS_REG_CONFIG1     0x01
#define ADS_REG_CONFIG2     0x02
#define ADS_REG_CONFIG3     0x03
#define ADS_REG_LOFF        0x04
#define ADS_REG_CH1SET      0x05
#define ADS_REG_CH2SET      0x06
#define ADS_REG_CH3SET      0x07
#define ADS_REG_CH4SET      0x08
#define ADS_REG_CH5SET      0x09
#define ADS_REG_CH6SET      0x0A
#define ADS_REG_CH7SET      0x0B
#define ADS_REG_CH8SET      0x0C
#define ADS_REG_RLD_SENSP   0x0D
#define ADS_REG_RLD_SENSN   0x0E
#define ADS_REG_LOFF_SENSP  0x0F
#define ADS_REG_LOFF_SENSN  0x10
#define ADS_REG_LOFF_FLIP   0x11
#define ADS_REG_GPIO        0x14
#define ADS_REG_PACE        0x15
#define ADS_REG_RESP        0x16
#define ADS_REG_CONFIG4     0x17
#define ADS_REG_WCT1        0x18
#define ADS_REG_WCT2        0x19

// ======================== CONSTANTS ==============================
#define NUM_CH           8
#define BAUD_RATE   250000
// SPI: ADS1298 uses Mode 1 (CPOL=0 CPHA=1), 4 MHz is safe
static const SPISettings adsSPI(4000000, MSBFIRST, SPI_MODE1);

// ======================== GAIN LOOKUP TABLE =======================
// ADS1298 CHnSET register bits [6:4] = PGA gain setting
// Index by PGA bits value → actual gain multiplier
static const uint8_t GAIN_TO_PGA_BITS[] = {
//  gain: 1     2     3     4     6     8     12
    0x10, 0x20, 0x30, 0x40, 0x00, 0x50, 0x60
};
static const uint8_t VALID_GAINS[] = { 1, 2, 3, 4, 6, 8, 12 };
#define NUM_VALID_GAINS  7
static volatile uint8_t currentGain = 6;  // default gain

// ======================== BINARY PACKET PROTOCOL =================
// Frame: [0xA5][0x5A][count][sample0][sample1]...[sampleN-1]
// Each sample: [4B timestamp_us LE][27B raw ADS1298 SPI data] = 31B
#define SYNC_BYTE_1        0xA5
#define SYNC_BYTE_2        0x5A
#define RAW_SPI_BYTES      27     // 3B status + 8×3B channels
#define SAMPLE_SIZE        (4 + RAW_SPI_BYTES)  // 31 bytes per sample
#define SAMPLES_PER_FRAME  10     // samples buffered before TCP send
#define FRAME_HEADER_SIZE  3      // 2B sync + 1B count
#define FRAME_MAX_SIZE     (FRAME_HEADER_SIZE + SAMPLES_PER_FRAME * SAMPLE_SIZE)

// ======================== WIFI / TCP =============================
#define WIFI_SSID       "ESP32_ECG"
#define WIFI_PASS       "12345678"     // min 8 chars for WPA2
#define TCP_PORT_DATA   3333           // ESP → PC  (ECG stream)
#define TCP_PORT_CMD    3334           // PC → ESP  (reserved)

static WiFiServer serverData(TCP_PORT_DATA);
static WiFiServer serverCmd(TCP_PORT_CMD);
static WiFiClient clientData;          // single data client
static WiFiClient clientCmd;           // single command client

// ======================== SLEEP CONFIG ===========================
#define SLEEP_TIMEOUT_MS  (2UL * 60UL * 1000UL)   // 2 minutes
static unsigned long bootTime            = 0;
static bool          clientEverConnected = false;
static unsigned long lastCountdownPrint  = 0;

// ======================== RING BUFFER ============================
// Decouples high-priority SPI reads from lower-priority TCP sends.
// Size must be power of 2 for fast masking.
#define RING_BUF_SAMPLES  64
#define RING_BUF_MASK     (RING_BUF_SAMPLES - 1)

typedef struct {
    uint32_t timestamp_us;
    uint8_t  raw[RAW_SPI_BYTES];
} __attribute__((packed)) Sample_t;

static volatile Sample_t ringBuffer[RING_BUF_SAMPLES];
static volatile uint16_t ringHead = 0;       // next write position (SPI task)
static volatile uint16_t ringTail = 0;       // next read  position (TCP task)
static volatile uint32_t ringOverflows = 0;  // debug: lost samples
static volatile uint32_t totalSamplesRead = 0;

// SPI transfer buffer (WORD_ALIGNED for DMA safety)
static uint8_t WORD_ALIGNED_ATTR spiTxBuf[RAW_SPI_BYTES];  // all zeros

// ======================== DRDY ISR & SEMAPHORE ===================
static volatile uint32_t isrTimestamp = 0;
static SemaphoreHandle_t drdySemaphore = NULL;

void IRAM_ATTR onDRDY() {
    isrTimestamp = (uint32_t)micros();   // capture timestamp at exact DRDY edge
    BaseType_t xHigherPriorityTaskWoken = pdFALSE;
    xSemaphoreGiveFromISR(drdySemaphore, &xHigherPriorityTaskWoken);
    if (xHigherPriorityTaskWoken) {
        portYIELD_FROM_ISR();
    }
}

// Frame buffer used by TCP send task
static uint8_t tcpFrameBuffer[FRAME_MAX_SIZE];

// ======================== ADS1298 LOW-LEVEL SPI ==================
static inline void adsCS(bool sel) {
    digitalWrite(PIN_CS_ADS, sel ? LOW : HIGH);
}

static void adsSendCmd(uint8_t cmd) {
    adsCS(true);
    SPI.beginTransaction(adsSPI);
    SPI.transfer(cmd);
    SPI.endTransaction();
    adsCS(false);
    delayMicroseconds(4);          // tCMD ≈ 4×tCLK
}

static uint8_t adsRReg(uint8_t reg) {
    adsCS(true);
    SPI.beginTransaction(adsSPI);
    SPI.transfer(0x20 | (reg & 0x1F));   // RREG opcode
    SPI.transfer(0x00);                   // read 1 register
    delayMicroseconds(5);                 // tSDECODE
    uint8_t v = SPI.transfer(0x00);
    SPI.endTransaction();
    adsCS(false);
    delayMicroseconds(2);
    return v;
}

static void adsWReg(uint8_t reg, uint8_t val) {
    adsCS(true);
    SPI.beginTransaction(adsSPI);
    SPI.transfer(0x40 | (reg & 0x1F));   // WREG opcode
    SPI.transfer(0x00);                   // write 1 register
    SPI.transfer(val);
    SPI.endTransaction();
    adsCS(false);
    delayMicroseconds(2);
}

// ======================== ADS1298 INIT ===========================
bool adsInit() {

    // ---------- 1. Power-up sequence ----------
    digitalWrite(PIN_ADS_START, LOW);
    digitalWrite(PIN_ADS_PWDN,  HIGH);     // chip active
    digitalWrite(PIN_ADS_RESET, HIGH);
    adsCS(false);                           // CS high
    delay(500);            // wait tPOR ≈ 2^18 / fCLK ≈ 128 ms @ 2.048 MHz

    // ---------- 2. Hardware reset pulse ----------
    digitalWrite(PIN_ADS_RESET, LOW);
    delayMicroseconds(10);                 // min 2×tCLK ≈ 1 µs
    digitalWrite(PIN_ADS_RESET, HIGH);
    delay(100);                            // reset recovery

    // ---------- 3. Stop continuous read (must before reg access) -
    adsSendCmd(ADS_CMD_SDATAC);
    delay(10);

    // ---------- 4. Verify device ID ----------
    uint8_t id = adsRReg(ADS_REG_ID);
    Serial.printf("ADS1298 ID: 0x%02X\n", id);
    // ADS1298 DEV_ID bits[7:5] should be in range 0x80..0xBF typically
    if (id == 0x00 || id == 0xFF) {
        Serial.println("ERROR: ADS1298 not responding!");
        return false;
    }

    // ---------- 5. Configure registers ----------

    // CONFIG1  (0x01)
    //   [7]   HR       = 1  (high-resolution mode)
    //   [6]   DAISY_EN = 0  (multiple readback; DAISY_IN=GND)
    //   [5]   CLK_EN   = 0  (osc clock output off)
    //   [2:0] DR       = 110 → 500 SPS in HR mode
    adsWReg(ADS_REG_CONFIG1, 0x86);

    // CONFIG2  (0x02) – no test signals
    adsWReg(ADS_REG_CONFIG2, 0x00);

    // CONFIG3  (0x03)
    //   [7]   PD_REFBUF  = 1  (internal reference ON)
    //   [6]   reserved   = 1  (must be 1)
    //   [5]   VREF_4V    = 0  (VREFP = 2.4 V)
    //   [3]   RLDREF_INT = 1  (RLD reference = internal AVDD/2)
    //   [2]   PD_RLD     = 1  (RLD amplifier ON → drives RLDOUT)
    adsWReg(ADS_REG_CONFIG3, 0xCC);
    delay(150);                            // wait for internal ref to settle

    // CONFIG4  – default
    adsWReg(ADS_REG_CONFIG4, 0x00);

    // LOFF – lead-off detection off
    adsWReg(ADS_REG_LOFF, 0x00);

    // Channel settings: all ON, Gain = 6, Normal electrode input
    //   [7]   PD       = 0     (channel ON)
    //   [6:4] PGA[2:0] = 000   (gain 6)
    //   [3]   reserved = 0
    //   [2:0] MUX[2:0] = 000   (normal electrode input)
    const uint8_t chCfg = 0x00;   // PGA=000→gain6, MUX=000→normal
    for (uint8_t r = ADS_REG_CH1SET; r <= ADS_REG_CH8SET; r++)
        adsWReg(r, chCfg);

    // RLD – driven from RA, LA, LL for common-mode rejection
    //   RLD_SENSP: CH2P (LA) + CH3P (LL) → bits [2:1] = 0x06
    //   RLD_SENSN: CH2N (RA)              → bit  [1]   = 0x02
    adsWReg(ADS_REG_RLD_SENSP, 0x06);
    adsWReg(ADS_REG_RLD_SENSN, 0x02);

    // Lead-off detection – disabled
    adsWReg(ADS_REG_LOFF_SENSP, 0x00);
    adsWReg(ADS_REG_LOFF_SENSN, 0x00);
    adsWReg(ADS_REG_LOFF_FLIP,  0x00);

    // ---------- 6. Wilson Central Terminal (WCT) ----------
    // WCT = (RA + LA + LL) / 3
    //
    // WCT1  (0x18)
    //   [7:4] aVF/aVL/aVR bits = 0000  (not used)
    //   [3]   PD_WCTA = 1  (WCTA ON)
    //   [2:0] WCTA    = 011 → CH2N (IN2N = RA)
    adsWReg(ADS_REG_WCT1, 0x0B);

    // WCT2  (0x19)
    //   [7]   PD_WCTC = 1
    //   [6:4] WCTC    = 100 → CH3P (IN3P = LL)
    //   [3]   PD_WCTB = 1
    //   [2:0] WCTB    = 010 → CH2P (IN2P = LA)
    adsWReg(ADS_REG_WCT2, 0xCA);

    // ---------- 7. Read-back key registers for debug ----------
    Serial.printf("CONFIG1 = 0x%02X\n", adsRReg(ADS_REG_CONFIG1));
    Serial.printf("CONFIG2 = 0x%02X\n", adsRReg(ADS_REG_CONFIG2));
    Serial.printf("CONFIG3 = 0x%02X\n", adsRReg(ADS_REG_CONFIG3));
    Serial.printf("CH1SET  = 0x%02X\n", adsRReg(ADS_REG_CH1SET));
    Serial.printf("WCT1    = 0x%02X\n", adsRReg(ADS_REG_WCT1));
    Serial.printf("WCT2    = 0x%02X\n", adsRReg(ADS_REG_WCT2));
    currentGain = 6;  // initial gain

    // ---------- 8. Start conversions ----------
    digitalWrite(PIN_ADS_START, HIGH);
    delay(10);

    // ---------- 9. RDATA mode (no RDATAC) ----------
    // Device stays in SDATAC state (default after reset / step 3).
    // Each sample is fetched explicitly via RDATA command in spiReadTask.
    // This avoids continuous-read desynchronisation and allows register
    // access at any time without SDATAC/RDATAC transitions.

    // ---------- 10. Attach DRDY interrupt (falling edge) ----------
    // NOTE: interrupt attached later in setup(), after tasks are created
    Serial.println("ADS1298 init OK – RDATA mode, 500 SPS streaming");
    return true;
}

// ======================== SET PGA GAIN ===========================
// Stops conversions, writes new gain to all CHnSET registers, restarts.
// In RDATA mode, register access is allowed without SDATAC/RDATAC
// transitions – we only need to stop and restart conversions.
// Returns true on success.
static bool adsSetGain(uint8_t gain) {
    // Find PGA bits for requested gain
    int idx = -1;
    for (int i = 0; i < NUM_VALID_GAINS; i++) {
        if (VALID_GAINS[i] == gain) { idx = i; break; }
    }
    if (idx < 0) return false;

    uint8_t pgaBits = GAIN_TO_PGA_BITS[idx];

    // Detach DRDY interrupt to prevent SPI conflicts
    detachInterrupt(digitalPinToInterrupt(PIN_ADS_DRDY));
    vTaskDelay(pdMS_TO_TICKS(5));  // let any in-progress SPI read finish

    // Stop conversions
    digitalWrite(PIN_ADS_START, LOW);
    delay(10);

    // In RDATA mode, device is not in RDATAC state, so register
    // access is allowed directly – no SDATAC needed.

    // Write new gain to all 8 channel registers
    // CHnSET: [7]=PD(0=ON), [6:4]=PGA, [3]=0, [2:0]=MUX(000=normal)
    uint8_t chCfg = pgaBits;  // PD=0(ON), MUX=000(normal)
    for (uint8_t r = ADS_REG_CH1SET; r <= ADS_REG_CH8SET; r++) {
        adsWReg(r, chCfg);
    }

    // Verify write
    uint8_t readBack = adsRReg(ADS_REG_CH1SET);
    Serial.printf("GAIN set: %d → CHnSET=0x%02X (readback=0x%02X)\n",
                  gain, chCfg, readBack);

    // Restart conversions
    digitalWrite(PIN_ADS_START, HIGH);
    delay(10);

    // Re-attach DRDY interrupt
    attachInterrupt(digitalPinToInterrupt(PIN_ADS_DRDY), onDRDY, FALLING);

    currentGain = gain;
    return (readBack == chCfg);
}

// ======================== ATOMIC SPI RDATA + READ 27 BYTES =======
// Sends RDATA opcode, waits tCMD, then block-transfers 27 data bytes.
// Entire sequence runs with interrupts disabled to prevent WiFi IRQ
// from inserting delays between bytes (spike elimination).
static void adsReadRdataAtomic(uint8_t* buf27) {
    portDISABLE_INTERRUPTS();

    digitalWrite(PIN_CS_ADS, LOW);
    SPI.beginTransaction(adsSPI);

    // 1) Send RDATA command (0x12)
    SPI.transfer(ADS_CMD_RDATA);

    // 2) Wait tCMD ≈ 4 × tCLK.  ADS1298 internal clock ≈ 2.048 MHz
    //    → 4 / 2.048 MHz ≈ 2 µs.  Use 4 µs for margin.
    delayMicroseconds(4);

    // 3) Clock out 27 bytes (3B status + 8×3B channel data)
    SPI.transferBytes(spiTxBuf, buf27, RAW_SPI_BYTES);

    SPI.endTransaction();
    digitalWrite(PIN_CS_ADS, HIGH);

    portENABLE_INTERRUPTS();
}

// ======================== SPI READ TASK ==========================
// High-priority task: waits for DRDY semaphore, sends RDATA command,
// reads ADS1298 data, validates status byte, stores in ring buffer.
static void spiReadTask(void* param) {
    (void)param;
    memset(spiTxBuf, 0, sizeof(spiTxBuf));
    Serial.println("SPI Read Task started (priority 5) – RDATA mode");

    while (true) {
        // Wait for DRDY signal from ISR (timeout 10 ms = 5× period @ 500 Hz)
        if (xSemaphoreTake(drdySemaphore, pdMS_TO_TICKS(10)) != pdTRUE) {
            continue;  // timeout – no DRDY, loop back
        }

        // Grab timestamp captured at exact DRDY falling edge
        uint32_t ts = isrTimestamp;

        // Send RDATA + atomic SPI read of all 27 bytes
        uint8_t rawBuf[RAW_SPI_BYTES];
        adsReadRdataAtomic(rawBuf);

        // --- Validate ADS1298 status byte ---
        // Per datasheet: status word bits [23:20] = 1100 (0xC0 mask on first byte)
        // If this doesn't match, the SPI read is misaligned → discard sample
        uint8_t statusHigh = rawBuf[0];
        if ((statusHigh & 0xF0) != 0xC0) {
            static uint32_t badStatusCount = 0;
            badStatusCount++;
            if ((badStatusCount & 0xFF) == 1) {  // print every 256th
                Serial.printf("WARN: Bad ADS status 0x%02X (total=%lu)\n",
                              statusHigh, badStatusCount);
            }
            continue;  // skip corrupt sample
        }

        // --- Store in ring buffer ---
        uint16_t nextHead = (ringHead + 1) & RING_BUF_MASK;
        if (nextHead == ringTail) {
            // Ring buffer full – overrun, drop oldest sample
            ringOverflows = (uint32_t)(ringOverflows + 1);
            ringTail = (ringTail + 1) & RING_BUF_MASK;
        }

        ringBuffer[ringHead].timestamp_us = ts;
        memcpy((void*)ringBuffer[ringHead].raw, rawBuf, RAW_SPI_BYTES);

        // Compiler memory barrier (ensures stores are visible to TCP task)
        __asm__ __volatile__("" ::: "memory");

        ringHead = nextHead;
        totalSamplesRead = (uint32_t)(totalSamplesRead + 1);
    }
}

// ======================== TCP SEND TASK ==========================
// Lower-priority task: drains ring buffer into TCP frame buffer,
// sends frame every SAMPLES_PER_FRAME samples.
static void tcpSendTask(void* param) {
    (void)param;
    Serial.println("TCP Send Task started (priority 2)");

    uint8_t frameSampleCount = 0;
    uint32_t lastStatsPrint = 0;

    while (true) {
        // --- Accept new TCP data client ---
        if (!clientData || !clientData.connected()) {
            WiFiClient newC = serverData.accept();
            if (newC) {
                clientData = newC;
                clientData.setNoDelay(true);
                clientEverConnected = true;
                frameSampleCount = 0;   // clean start for new client
                Serial.printf("Data client connected: %s\n",
                              clientData.remoteIP().toString().c_str());
            }
        }

        // --- Accept new TCP cmd client ---
        if (!clientCmd || !clientCmd.connected()) {
            WiFiClient newC = serverCmd.accept();
            if (newC) {
                clientCmd = newC;
                clientCmd.setNoDelay(true);
                Serial.printf("Cmd  client connected: %s\n",
                              clientCmd.remoteIP().toString().c_str());
            }
        }

        // --- Handle commands from cmd port ---
        if (clientCmd && clientCmd.connected() && clientCmd.available()) {
            static char cmdBuf[64];
            static uint8_t cmdLen = 0;
            while (clientCmd.available() && cmdLen < sizeof(cmdBuf) - 1) {
                char c = clientCmd.read();
                if (c == '\n' || c == '\r') {
                    if (cmdLen > 0) {
                        cmdBuf[cmdLen] = '\0';
                        String cmd(cmdBuf);
                        Serial.printf("CMD received: %s\n", cmdBuf);

                        if (cmd.startsWith("GAIN:")) {
                            int reqGain = cmd.substring(5).toInt();
                            if (adsSetGain((uint8_t)reqGain)) {
                                char resp[64];
                                snprintf(resp, sizeof(resp), "GAIN_OK:%d\n", currentGain);
                                clientCmd.print(resp);
                                Serial.printf("GAIN_OK: %d\n", currentGain);
                            } else {
                                char resp[64];
                                snprintf(resp, sizeof(resp), "GAIN_ERR:invalid gain %d\n", reqGain);
                                clientCmd.print(resp);
                                Serial.printf("GAIN_ERR: invalid gain %d\n", reqGain);
                            }
                        } else {
                            clientCmd.print("ERR:unknown command\n");
                        }
                        cmdLen = 0;
                    }
                } else {
                    cmdBuf[cmdLen++] = c;
                }
            }
            if (cmdLen >= sizeof(cmdBuf) - 1) cmdLen = 0;  // overflow reset
        }

        // --- Collect samples from ring buffer into frame ---
        bool gotSample = false;
        while (ringTail != ringHead && frameSampleCount < SAMPLES_PER_FRAME) {
            uint16_t tail = ringTail;

            uint32_t ts = ringBuffer[tail].timestamp_us;
            uint8_t* dst = &tcpFrameBuffer[FRAME_HEADER_SIZE +
                                           frameSampleCount * SAMPLE_SIZE];

            // Timestamp (little-endian 4 bytes)
            dst[0] = (uint8_t)(ts         & 0xFF);
            dst[1] = (uint8_t)((ts >> 8)  & 0xFF);
            dst[2] = (uint8_t)((ts >> 16) & 0xFF);
            dst[3] = (uint8_t)((ts >> 24) & 0xFF);

            // Copy raw SPI data (27 bytes: 3B status + 8×3B channels)
            memcpy(dst + 4, (const void*)ringBuffer[tail].raw, RAW_SPI_BYTES);

            // Memory barrier before advancing tail
            __asm__ __volatile__("" ::: "memory");
            ringTail = (tail + 1) & RING_BUF_MASK;

            frameSampleCount++;
            gotSample = true;
        }

        // --- Send frame when full ---
        if (frameSampleCount >= SAMPLES_PER_FRAME) {
            tcpFrameBuffer[0] = SYNC_BYTE_1;
            tcpFrameBuffer[1] = SYNC_BYTE_2;
            tcpFrameBuffer[2] = frameSampleCount;

            if (clientData && clientData.connected()) {
                size_t frameSize = FRAME_HEADER_SIZE +
                                   frameSampleCount * SAMPLE_SIZE;
                size_t written = clientData.write(tcpFrameBuffer, frameSize);
                if (written != frameSize) {
                    Serial.printf("WARN: TCP wrote %d/%d bytes\n",
                                  written, frameSize);
                }
            }
            frameSampleCount = 0;
        }

        // --- Debug stats every 5 seconds ---
        if (millis() - lastStatsPrint >= 5000) {
            lastStatsPrint = millis();
            uint16_t used = (ringHead - ringTail) & RING_BUF_MASK;
            Serial.printf("[STATS] samples=%lu ring=%d/%d overflows=%lu\n",
                          totalSamplesRead, used,
                          RING_BUF_SAMPLES, ringOverflows);
        }

        // If no data was available, yield briefly to avoid busy-spin
        if (!gotSample) {
            vTaskDelay(pdMS_TO_TICKS(1));
        }
    }
}

// ======================== WIFI INIT ==============================
static void wifiInit() {
    WiFi.mode(WIFI_AP);
    WiFi.softAP(WIFI_SSID, WIFI_PASS);
    delay(200);   // let DHCP settle

    IPAddress ip = WiFi.softAPIP();
    Serial.printf("WiFi AP \"%s\" started – IP: %s\n", WIFI_SSID, ip.toString().c_str());

    serverData.begin();
    serverData.setNoDelay(true);
    Serial.printf("TCP data   server on port %d\n", TCP_PORT_DATA);

    serverCmd.begin();
    serverCmd.setNoDelay(true);
    Serial.printf("TCP cmd    server on port %d (reserved)\n", TCP_PORT_CMD);
}

// ======================== SETUP ==================================
void setup() {
    Serial.begin(BAUD_RATE);
    delay(4000);                           // wait for USB-CDC ready
    Serial.println("\n=== ESP32-C6  12-Lead ECG  v11 (RDATA mode) ===");

    // --- pin modes ---
    pinMode(PIN_CS_ADS,    OUTPUT);
    pinMode(PIN_CS_MPU,    OUTPUT);
    pinMode(PIN_ADS_DRDY,  INPUT);
    pinMode(PIN_ADS_PWDN,  OUTPUT);
    pinMode(PIN_ADS_RESET, OUTPUT);
    pinMode(PIN_ADS_START, OUTPUT);
    pinMode(PIN_MPU_INT,   INPUT);

    // disable MPU9250 – CS stays HIGH
    digitalWrite(PIN_CS_MPU, HIGH);
    // deselect ADS
    digitalWrite(PIN_CS_ADS, HIGH);

    // init SPI (CS managed manually → pass -1)
    SPI.begin(PIN_SPI_CLK, PIN_SPI_MISO, PIN_SPI_MOSI, -1);

    // record boot time for sleep-timeout logic
    bootTime = millis();

    // init WiFi AP + TCP servers
    wifiInit();

    Serial.printf("Sleep timeout: %lu s – waiting for client...\n",
                  SLEEP_TIMEOUT_MS / 1000);

    // init ADS1298 (DRDY interrupt NOT attached yet)
    if (!adsInit()) {
        Serial.println("FATAL: ADS1298 init failed – halting.");
        while (true) delay(1000);
    }

    // Create binary semaphore for DRDY signalling
    drdySemaphore = xSemaphoreCreateBinary();

    // Create SPI read task – HIGH priority (must not miss any DRDY)
    // ESP32-C6 is single-core, but priority preemption still works:
    // SPI task (prio 5) preempts TCP task (prio 2) whenever DRDY fires.
    xTaskCreate(
        spiReadTask,
        "spi_read",
        4096,       // stack size
        NULL,       // param
        5,          // priority: HIGH (above WiFi tasks ~1, above TCP task)
        NULL        // task handle (not needed)
    );

    // Create TCP send task – LOWER priority
    xTaskCreate(
        tcpSendTask,
        "tcp_send",
        8192,       // stack size (more for WiFi/TCP operations)
        NULL,
        2,          // priority: above idle, below SPI read
        NULL
    );

    // NOW attach DRDY interrupt (tasks are ready to receive)
    attachInterrupt(digitalPinToInterrupt(PIN_ADS_DRDY), onDRDY, FALLING);

    Serial.println("All tasks started – system running.");
}

// ======================== LOOP ===================================
// On v11, loop() only handles the deep-sleep timeout logic.
// All ECG data flow is handled by FreeRTOS tasks (spiReadTask, tcpSendTask).
void loop() {

    // --- deep-sleep check: no client within timeout → sleep ----------
    if (!clientEverConnected) {
        unsigned long elapsed = millis() - bootTime;

        // print countdown every 10 seconds
        if (millis() - lastCountdownPrint >= 10000) {
            lastCountdownPrint = millis();
            unsigned long remaining_s = (elapsed < SLEEP_TIMEOUT_MS)
                                        ? (SLEEP_TIMEOUT_MS - elapsed) / 1000
                                        : 0;
            Serial.printf("No client yet – sleep in %lu s\n", remaining_s);
        }

        if (elapsed >= SLEEP_TIMEOUT_MS) {
            Serial.println("No client connected within timeout – entering deep sleep...");
            Serial.flush();
            // stop ADS1298 conversions & WiFi before sleep
            adsSendCmd(ADS_CMD_STOP);
            WiFi.softAPdisconnect(true);
            WiFi.mode(WIFI_OFF);
            delay(100);
            // enter deep sleep (only reset will wake the chip)
            esp_deep_sleep_start();
            // execution never reaches here
        }
    }

    // loop() is mostly idle now – real work is in FreeRTOS tasks
    vTaskDelay(pdMS_TO_TICKS(100));
}