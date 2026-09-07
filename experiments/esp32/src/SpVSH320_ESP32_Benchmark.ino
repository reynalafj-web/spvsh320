/*
  SpVSH-320-FR2 ESP32 Arduino Benchmark

  Purpose:
    - Measure implementation metrics for SpVSH-320-FR2 on ESP32 using Arduino IDE.
    - Reports timing, derived cycles, CPB, throughput, heap delta, and digest sample.
    - Supports optional 1-hour endurance test with 1024-byte payload for FNIRSI FNB58 measurement.

  Board target:
    - ESP32 DevKit / ESP32 Arduino Core
    - Recommended CPU clock: 240 MHz

  Notes:
    - Wi-Fi and Bluetooth are disabled in setup to reduce power variability.
    - Do not print inside the measured hashing loop.
    - FNIRSI FNB58 measures board-level USB input power, not isolated cryptographic-core power.

  Algorithm:
    - SpVSH-320-FR2, 320-bit state, 256-bit digest.
    - Two output-finalization rounds before first output word.
*/
#pragma GCC optimize ("O3,unroll-loops")
#include <Arduino.h>
#include <WiFi.h>
#include "esp_wifi.h"
#include "esp_bt.h"

// ============================================================
// User configuration
// ============================================================
static const bool RUN_SIZE_BENCHMARK = true;
static const bool RUN_ENDURANCE_TEST = false;       // Set true for FNIRSI FNB58 1-hour test.
static const uint32_t SERIAL_BAUD = 115200;

static const uint32_t QUICK_ITERATIONS = 1000;      // Per message size. Increase if timing is too noisy.
static const uint32_t WARMUP_ITERATIONS = 100;

static const size_t ENDURANCE_PAYLOAD_BYTES = 1024;
static const uint32_t ENDURANCE_SECONDS = 3600;     // 1 hour.
static const uint32_t ENDURANCE_PROGRESS_SECONDS = 60;

static const uint32_t FINALIZATION_ROUNDS = 2;      // SpVSH-320-FR2.

// ============================================================
// SpVSH-320-FR2 core implementation
// ============================================================
struct PentaState {
  uint64_t s[5];
};

static constexpr uint64_t P1 = 0xFFFFFFFFFFFFFFC5ULL; // 2^64 - 59
static constexpr uint64_t C1 = 59;
static constexpr uint64_t P2 = 0xFFFFFFFFFFFFFFADULL; // 2^64 - 83
static constexpr uint64_t C2 = 83;

static const uint16_t SMALL_PRIMES[64] = {
    3,   5,   7,  11,  13,  17,  19,  23,  29,  31,  37,  41,  43,  47,  53,  59,
   61,  67,  71,  73,  79,  83,  89,  97, 101, 103, 107, 109, 113, 127, 131, 137,
  139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227,
  229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313
};

static inline uint64_t rotr64(uint64_t x, uint8_t r) {
  return (x >> r) | (x << (64 - r));
}

static inline void mul64x64_128(uint64_t a, uint64_t b, uint64_t *res_hi, uint64_t *res_lo) {
  uint32_t a_lo = static_cast<uint32_t>(a);
  uint32_t a_hi = static_cast<uint32_t>(a >> 32);
  uint32_t b_lo = static_cast<uint32_t>(b);
  uint32_t b_hi = static_cast<uint32_t>(b >> 32);

  uint64_t t = static_cast<uint64_t>(a_lo) * b_lo;
  uint32_t w3 = static_cast<uint32_t>(t);
  uint32_t k = static_cast<uint32_t>(t >> 32);

  t = static_cast<uint64_t>(a_hi) * b_lo + k;
  uint32_t w2 = static_cast<uint32_t>(t);
  uint32_t w1 = static_cast<uint32_t>(t >> 32);

  t = static_cast<uint64_t>(a_lo) * b_hi + w2;
  k = static_cast<uint32_t>(t >> 32);

  *res_lo = (static_cast<uint64_t>(static_cast<uint32_t>(t)) << 32) | w3;
  *res_hi = static_cast<uint64_t>(a_hi) * b_hi + w1 + k;
}

static inline uint64_t MulModP(uint64_t a, uint64_t b, uint64_t P, uint64_t C) {
  uint64_t hi, lo;
  mul64x64_128(a, b, &hi, &lo);

  uint64_t h_lo = static_cast<uint32_t>(hi);
  uint64_t h_hi = (hi >> 32);

  uint64_t p1 = h_lo * C;
  uint64_t p2 = h_hi * C;

  uint64_t l2 = p1 + (p2 << 32);
  uint64_t h2 = (p2 >> 32) + (l2 < p1 ? 1 : 0);

  uint64_t sum_lo = l2 + lo;
  uint64_t carry = (sum_lo < lo) ? 1 : 0;
  uint64_t sum_hi = h2 + carry;

  uint64_t final_res = sum_lo + sum_hi * C;

  if (final_res < sum_lo) final_res += C;
  if (final_res >= P) final_res -= P;

  return final_res;
}

static inline uint64_t MulModP1(uint64_t a, uint64_t b) {
  return MulModP(a, b, P1, C1);
}

static inline uint64_t MulModP2(uint64_t a, uint64_t b) {
  return MulModP(a, b, P2, C2);
}

static inline uint64_t GET_U64_LE(const uint8_t* ptr) {
  uint64_t val;
  memcpy(&val, ptr, 8);
  return val;
}

static inline void PUT_U64_LE(uint64_t x, uint8_t* out) {
  memcpy(out, &x, 8);
}

static inline void PermuteSpVSH(PentaState &state) {
  uint64_t trigger = state.s[0] ^ state.s[2] ^ state.s[4] ^ 0x9E3779B97F4A7C15ULL;

  uint64_t prod_p1 = 1;
  uint64_t prod_p2 = 1;

  while (trigger != 0) {
    int idx = __builtin_ctzll(trigger);
    uint64_t p = SMALL_PRIMES[idx];

    prod_p1 = MulModP1(prod_p1, p);
    prod_p2 = MulModP2(prod_p2, p);

    trigger &= trigger - 1;
  }

  state.s[0] = MulModP1(MulModP1(state.s[0], state.s[0]), prod_p1);
  state.s[1] = MulModP2(MulModP2(state.s[1], state.s[1]), prod_p2);
  state.s[2] = MulModP1(MulModP1(state.s[2], state.s[2]), prod_p1);
  state.s[3] = MulModP2(MulModP2(state.s[3], state.s[3]), prod_p2);
  state.s[4] = MulModP1(MulModP1(state.s[4], state.s[4]), prod_p1);

  uint64_t t0 = state.s[0];
  uint64_t t1 = state.s[1];
  uint64_t t2 = state.s[2];
  uint64_t t3 = state.s[3];
  uint64_t t4 = state.s[4];

  state.s[0] = t0 ^ rotr64(t1, 17);
  state.s[1] = t1 ^ rotr64(t2, 11);
  state.s[2] = t2 ^ rotr64(t3, 19);
  state.s[3] = t3 ^ rotr64(t4, 13);
  state.s[4] = t4 ^ rotr64(t0, 23);
}

static inline PentaState InitSpVSHState() {
  PentaState state;
  uint64_t seed64 = 0x0123456789ABCDEFULL;

  state.s[0] = seed64 + 0xDEADBEEFCAFEBABEULL;
  state.s[1] = seed64 ^ 0x13198A2E03707344ULL;
  state.s[2] = seed64 + 0x9E3779B97F4A7C15ULL;
  state.s[3] = seed64 ^ 0xBF58476D1CE4E5B9ULL;
  state.s[4] = seed64 + 0x94D049BB133111EBULL;

  for (int i = 0; i < 5; i++) {
    if (state.s[i] == 0) state.s[i] = 3;
  }

  PermuteSpVSH(state);
  return state;
}

static inline void ApplyFinalizationRounds(PentaState& state, uint32_t rounds) {
  for (uint32_t r = 0; r < rounds; ++r) {
    PermuteSpVSH(state);
  }
}

void SpVSH320_FR2(const uint8_t *in, size_t len, uint8_t *out) {
  PentaState state = InitSpVSHState();

  size_t length = len;
  while (length >= 8) {
    state.s[0] ^= GET_U64_LE(in);
    PermuteSpVSH(state);
    in += 8;
    length -= 8;
  }

  uint8_t buffer[8] = {0};
  if (length > 0) {
    memcpy(buffer, in, length);
  }

  buffer[length] = 0x80;
  state.s[0] ^= GET_U64_LE(buffer);
  PermuteSpVSH(state);

  ApplyFinalizationRounds(state, FINALIZATION_ROUNDS);

  uint64_t h[4];
  for (int i = 0; i < 4; ++i) {
    h[i] = state.s[0];
    if (i < 3) {
      PermuteSpVSH(state);
    }
  }

  PUT_U64_LE(h[0], out + 0);
  PUT_U64_LE(h[1], out + 8);
  PUT_U64_LE(h[2], out + 16);
  PUT_U64_LE(h[3], out + 24);
}

// ============================================================
// Benchmark utilities
// ============================================================
static uint8_t payload1024[1024];
static uint8_t digest[32];

static void fillPayload(uint8_t* buf, size_t len, uint32_t domain) {
  uint32_t x = 0x9E3779B9UL ^ domain;
  for (size_t i = 0; i < len; ++i) {
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    buf[i] = static_cast<uint8_t>(x & 0xFF);
  }
}

static void printDigestHex(const uint8_t* d, size_t len) {
  for (size_t i = 0; i < len; ++i) {
    if (d[i] < 16) Serial.print('0');
    Serial.print(d[i], HEX);
  }
}

static uint32_t chooseIterations(size_t msgLen) {
  if (msgLen <= 16) return QUICK_ITERATIONS * 2;
  if (msgLen <= 128) return QUICK_ITERATIONS;
  if (msgLen <= 256) return max<uint32_t>(250, QUICK_ITERATIONS / 2);
  return max<uint32_t>(100, QUICK_ITERATIONS / 4);
}

static void runOneSizeBenchmark(size_t msgLen) {
  uint8_t* msg = payload1024;
  fillPayload(msg, msgLen, static_cast<uint32_t>(msgLen));

  for (uint32_t i = 0; i < WARMUP_ITERATIONS; ++i) {
    SpVSH320_FR2(msg, msgLen, digest);
  }

  const uint32_t iterations = chooseIterations(msgLen);
  const uint32_t heapBefore = ESP.getFreeHeap();
  const uint64_t startUs = esp_timer_get_time();

  for (uint32_t i = 0; i < iterations; ++i) {
    msg[0] ^= static_cast<uint8_t>(i);  // Prevent overly aggressive compiler assumptions.
    SpVSH320_FR2(msg, msgLen, digest);
    msg[0] ^= static_cast<uint8_t>(i);
  }

  const uint64_t endUs = esp_timer_get_time();
  const uint32_t heapAfter = ESP.getFreeHeap();

  const double totalUs = static_cast<double>(endUs - startUs);
  const double avgUs = totalUs / static_cast<double>(iterations);
  const double cpuMHz = static_cast<double>(ESP.getCpuFreqMHz());
  const double cyclesPerHash = avgUs * cpuMHz;
  const double cpb = cyclesPerHash / static_cast<double>(msgLen == 0 ? 1 : msgLen);
  const double throughputKBs = (static_cast<double>(msgLen) * static_cast<double>(iterations) * 1000000.0) / (totalUs * 1024.0);
  const int32_t heapDelta = static_cast<int32_t>(heapBefore) - static_cast<int32_t>(heapAfter);

  Serial.print("SpVSH-320-FR2,ESP32,");
  Serial.print(ESP.getCpuFreqMHz());
  Serial.print(",");
  Serial.print(msgLen);
  Serial.print(",");
  Serial.print(iterations);
  Serial.print(",");
  Serial.print(avgUs, 3);
  Serial.print(",");
  Serial.print(cyclesPerHash, 2);
  Serial.print(",");
  Serial.print(cpb, 2);
  Serial.print(",");
  Serial.print(throughputKBs, 2);
  Serial.print(",");
  Serial.print(heapBefore);
  Serial.print(",");
  Serial.print(heapAfter);
  Serial.print(",");
  Serial.print(heapDelta);
  Serial.print(",");
  printDigestHex(digest, 8);  // Short digest prefix for log sanity.
  Serial.println();
}

static void runSizeBenchmark() {
  Serial.println();
  Serial.println("=== SpVSH-320-FR2 ESP32 Size Benchmark ===");
  Serial.println("algorithm,board,cpu_clock_mhz,msg_size_bytes,iterations,avg_time_us,cycles_per_hash,cycles_per_byte,throughput_kb_s,heap_before,heap_after,heap_delta,digest_prefix64_hex");

  const size_t sizes[] = {0, 1, 16, 64, 128, 242, 256, 512, 1024};
  for (size_t i = 0; i < sizeof(sizes) / sizeof(sizes[0]); ++i) {
    runOneSizeBenchmark(sizes[i]);
    delay(200);
  }

  Serial.println("=== Size benchmark complete ===");
}

static void runEnduranceTest() {
  Serial.println();
  Serial.println("=== SpVSH-320-FR2 1-Hour Endurance Test ===");
  Serial.println("Use FNIRSI FNB58 to record active board-level energy during this phase.");
  Serial.println("Avoid serial monitor disconnection during the run.");
  Serial.print("payload_size_bytes=");
  Serial.println(ENDURANCE_PAYLOAD_BYTES);
  Serial.print("target_duration_seconds=");
  Serial.println(ENDURANCE_SECONDS);
  Serial.println("Starting in 10 seconds. Record idle power before this point.");
  delay(10000);

  fillPayload(payload1024, ENDURANCE_PAYLOAD_BYTES, 0x320F0B58UL);

  uint64_t totalHashes = 0;
  const uint64_t startUs = esp_timer_get_time();
  const uint64_t durationUs = static_cast<uint64_t>(ENDURANCE_SECONDS) * 1000000ULL;
  uint64_t nextProgressUs = static_cast<uint64_t>(ENDURANCE_PROGRESS_SECONDS) * 1000000ULL;

  while ((esp_timer_get_time() - startUs) < durationUs) {
    // No serial printing inside the tight loop except periodic progress below.
    SpVSH320_FR2(payload1024, ENDURANCE_PAYLOAD_BYTES, digest);
    ++totalHashes;

    const uint64_t elapsedUs = esp_timer_get_time() - startUs;
    if (elapsedUs >= nextProgressUs) {
      Serial.print("progress_seconds=");
      Serial.print(elapsedUs / 1000000ULL);
      Serial.print(",total_hashes=");
      Serial.println((unsigned long long)totalHashes);
      nextProgressUs += static_cast<uint64_t>(ENDURANCE_PROGRESS_SECONDS) * 1000000ULL;
    }
  }

  const uint64_t endUs = esp_timer_get_time();
  const double elapsedSec = static_cast<double>(endUs - startUs) / 1000000.0;
  const double avgTimeUs = (elapsedSec * 1000000.0) / static_cast<double>(totalHashes);
  const double totalBytes = static_cast<double>(totalHashes) * static_cast<double>(ENDURANCE_PAYLOAD_BYTES);
  const double throughputKBs = totalBytes / elapsedSec / 1024.0;
  const double cyclesPerHash = avgTimeUs * static_cast<double>(ESP.getCpuFreqMHz());
  const double cpb = cyclesPerHash / static_cast<double>(ENDURANCE_PAYLOAD_BYTES);

  Serial.println();
  Serial.println("=== Endurance result CSV ===");
  Serial.println("algorithm,board,cpu_clock_mhz,payload_size_bytes,duration_seconds,total_hashes,total_bytes_processed,avg_time_us,cycles_per_hash,cycles_per_byte,throughput_kb_s,digest_prefix64_hex");
  Serial.print("SpVSH-320,ESP32,");
  Serial.print(ESP.getCpuFreqMHz());
  Serial.print(",");
  Serial.print(ENDURANCE_PAYLOAD_BYTES);
  Serial.print(",");
  Serial.print(elapsedSec, 3);
  Serial.print(",");
  Serial.print((unsigned long long)totalHashes);
  Serial.print(",");
  Serial.print((unsigned long long)(totalHashes * ENDURANCE_PAYLOAD_BYTES));
  Serial.print(",");
  Serial.print(avgTimeUs, 3);
  Serial.print(",");
  Serial.print(cyclesPerHash, 2);
  Serial.print(",");
  Serial.print(cpb, 2);
  Serial.print(",");
  Serial.print(throughputKBs, 2);
  Serial.print(",");
  printDigestHex(digest, 8);
  Serial.println();
  Serial.println("=== Endurance test complete ===");
}

void disableRadiosForMeasurement() {
  WiFi.mode(WIFI_OFF);
  esp_wifi_stop();
  btStop();
  esp_bt_controller_disable();
}

void setup() {
  Serial.begin(SERIAL_BAUD);
  delay(2000);

  disableRadiosForMeasurement();
  setCpuFrequencyMhz(240);
  delay(500);

  fillPayload(payload1024, sizeof(payload1024), 0xA5A5A5A5UL);

  Serial.println("SpVSH-320-FR2 ESP32 Arduino Benchmark");
  Serial.print("CPU MHz: ");
  Serial.println(ESP.getCpuFreqMHz());
  Serial.print("Free heap: ");
  Serial.println(ESP.getFreeHeap());
  Serial.print("Finalization rounds: ");
  Serial.println(FINALIZATION_ROUNDS);

  SpVSH320_FR2((const uint8_t*)"abc", 3, digest);
  Serial.print("KAT digest prefix for 'abc': ");
  printDigestHex(digest, 16);
  Serial.println();

  if (RUN_SIZE_BENCHMARK) {
    runSizeBenchmark();
  }

  if (RUN_ENDURANCE_TEST) {
    runEnduranceTest();
  }

  Serial.println("All requested benchmarks completed.");
}

void loop() {
  delay(1000);
}
