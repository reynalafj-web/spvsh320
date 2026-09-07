/*
  SpVSH-320-FR2 ESP32 Arduino Benchmark

  Purpose:
    - Measure implementation metrics for SpVSH-320-FR2 on ESP32 using Arduino IDE.
    - Reports timing, derived cycles, CPB, throughput, heap delta, and digest sample.
    - Runs a 1-hour endurance test with 1024-byte payload for FNIRSI FNB58 measurement.

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
#include "esp_timer.h"
#include <WiFi.h>
#include "esp_wifi.h"
#include "esp_bt.h"

// ============================================================
// User configuration
// ============================================================
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
// Standardized ESP32 Benchmark Harness
// ============================================================
// Metrics:
//   - Throughput: derived from processed bytes and elapsed time.
//   - Energy consumption: computed only when BENCH_AVG_POWER_MW is filled
//     using external FNIRSI FNB58 measurement. ESP32 cannot measure board
//     energy directly without an external sensor.
//   - Power consumption: BENCH_AVG_POWER_MW should be taken from FNIRSI FNB58.
//   - ROM usage: ESP.getSketchSize().
//   - RAM usage: ESP.getFreeHeap(), ESP.getMinFreeHeap(), ESP.getMaxAllocHeap().
//
// Recommended protocol:
//   1. Flash only one algorithm sketch at a time.
//   2. Use the same ESP32 board, USB cable, power source, Arduino core, and CPU clock.
//   3. Disable Serial Monitor during endurance measurement when using FNIRSI.
//   4. Record average power from FNIRSI and insert it into BENCH_AVG_POWER_MW.
//   5. Re-run to obtain energy_per_hash_mJ and total_energy_mJ in the CSV output.

#ifndef ESP32
#error "This benchmark is intended for ESP32 boards."
#endif

static const char* BENCH_ALGORITHM = "SpVSH-320";
static const uint32_t BENCH_CPU_MHZ = 240;
static const uint32_t BENCH_SERIAL_BAUD = 115200;

static const uint32_t BENCH_WARMUP_ITERATIONS = 5;
static const uint32_t BENCH_MAX_ITERATIONS = 1000;
static const uint32_t BENCH_MIN_DURATION_MS = 1000;

static const bool RUN_SIZE_BENCHMARK = true;
static const bool RUN_ENDURANCE_TEST = false;
static const uint32_t ENDURANCE_SECONDS = 3600;
static const size_t ENDURANCE_PAYLOAD_BYTES = 1024;
static const uint32_t ENDURANCE_PROGRESS_SECONDS = 300;

// Fill this value using the average power shown by FNIRSI FNB58.
// Example: if FNIRSI reports 421.7 mW, set BENCH_AVG_POWER_MW = 421.7.
// Leave as 0.0 when no external power measurement is available.
static const double BENCH_AVG_POWER_MW = 0.0;

static const size_t BENCH_DIGEST_BYTES = 32;
static const size_t BENCH_MAX_PAYLOAD_BYTES = 1024;

static uint8_t bench_payload[BENCH_MAX_PAYLOAD_BYTES];
static uint8_t bench_digest[BENCH_DIGEST_BYTES];

static const size_t BENCH_SIZES[] = {16, 64, 128, 242, 256, 512, 1024};

static inline uint64_t benchNowUs() {
  return (uint64_t)esp_timer_get_time();
}

static void fillPayload(uint8_t* buf, size_t len, uint32_t seed) {
  uint32_t x = seed;
  for (size_t i = 0; i < len; ++i) {
    x = 1664525UL * x + 1013904223UL;
    buf[i] = static_cast<uint8_t>((x >> 24) ^ (i * 31U));
  }
}

static void printDigestPrefixHex(const uint8_t* d, size_t n) {
  static const char HEX_CHARS[] = "0123456789abcdef";
  for (size_t i = 0; i < n; ++i) {
    Serial.print(HEX_CHARS[(d[i] >> 4) & 0x0F]);
    Serial.print(HEX_CHARS[d[i] & 0x0F]);
  }
}

static inline void runHashFunction(const uint8_t* msg, size_t len, uint8_t out[32]) {
  SpVSH320_FR2(msg, len, out);
}

static void disableRadiosForMeasurement() {
  WiFi.mode(WIFI_OFF);
  WiFi.disconnect(true);
  btStop();
}

static void printSystemMetadata() {
  Serial.println();
  Serial.println("============================================================");
  Serial.print(BENCH_ALGORITHM);
  Serial.println(" Standardized ESP32 Benchmark");
  Serial.println("============================================================");
  Serial.print("cpu_mhz=");
  Serial.println(ESP.getCpuFreqMHz());
  Serial.print("sketch_size_bytes=");
  Serial.println(ESP.getSketchSize());
  Serial.print("free_sketch_space_bytes=");
  Serial.println(ESP.getFreeSketchSpace());
  Serial.print("flash_chip_size_bytes=");
  Serial.println(ESP.getFlashChipSize());
  Serial.print("free_heap_start_bytes=");
  Serial.println(ESP.getFreeHeap());
  Serial.print("min_free_heap_start_bytes=");
  Serial.println(ESP.getMinFreeHeap());
  Serial.print("max_alloc_heap_start_bytes=");
  Serial.println(ESP.getMaxAllocHeap());
  Serial.print("configured_external_avg_power_mw=");
  Serial.println(BENCH_AVG_POWER_MW, 3);
  Serial.println();
}

static void runOneSizeBenchmark(size_t msgSize) {
  fillPayload(bench_payload, msgSize, 0xA5A5A5A5UL + (uint32_t)msgSize);

  for (uint32_t i = 0; i < BENCH_WARMUP_ITERATIONS; ++i) {
    bench_payload[0] ^= static_cast<uint8_t>(i);
    runHashFunction(bench_payload, msgSize, bench_digest);
    bench_payload[0] ^= static_cast<uint8_t>(i);
  }

  const uint32_t heapBefore = ESP.getFreeHeap();
  const uint32_t minHeapBefore = ESP.getMinFreeHeap();
  const uint32_t maxAllocBefore = ESP.getMaxAllocHeap();

  const uint64_t startUs = benchNowUs();
  uint64_t elapsedUs = 0;
  uint32_t iterations = 0;

  do {
    bench_payload[0] ^= static_cast<uint8_t>(iterations & 0xFFU);
    runHashFunction(bench_payload, msgSize, bench_digest);
    bench_payload[0] ^= static_cast<uint8_t>(iterations & 0xFFU);
    ++iterations;
    elapsedUs = benchNowUs() - startUs;
  } while (iterations < BENCH_MAX_ITERATIONS &&
           elapsedUs < ((uint64_t)BENCH_MIN_DURATION_MS * 1000ULL));

  const uint32_t heapAfter = ESP.getFreeHeap();
  const uint32_t minHeapAfter = ESP.getMinFreeHeap();
  const uint32_t maxAllocAfter = ESP.getMaxAllocHeap();

  const double elapsedSec = static_cast<double>(elapsedUs) / 1000000.0;
  const double avgTimeUs = static_cast<double>(elapsedUs) / static_cast<double>(iterations);
  const double cyclesPerHash = avgTimeUs * static_cast<double>(ESP.getCpuFreqMHz());
  const double cyclesPerByte = cyclesPerHash / static_cast<double>(msgSize);
  const double totalBytes = static_cast<double>(iterations) * static_cast<double>(msgSize);
  const double throughputKBs = totalBytes / elapsedSec / 1024.0;

  const double totalEnergyMJ = (BENCH_AVG_POWER_MW > 0.0) ? (BENCH_AVG_POWER_MW * elapsedSec) : 0.0;
  const double energyPerHashMJ = (BENCH_AVG_POWER_MW > 0.0) ? (totalEnergyMJ / static_cast<double>(iterations)) : 0.0;
  const double energyPerByteUJ = (BENCH_AVG_POWER_MW > 0.0) ? ((totalEnergyMJ * 1000.0) / totalBytes) : 0.0;

  Serial.print(BENCH_ALGORITHM);
  Serial.print(",ESP32,");
  Serial.print(ESP.getCpuFreqMHz());
  Serial.print(",");
  Serial.print(msgSize);
  Serial.print(",");
  Serial.print(iterations);
  Serial.print(",");
  Serial.print(elapsedSec, 6);
  Serial.print(",");
  Serial.print(avgTimeUs, 3);
  Serial.print(",");
  Serial.print(cyclesPerHash, 3);
  Serial.print(",");
  Serial.print(cyclesPerByte, 3);
  Serial.print(",");
  Serial.print(throughputKBs, 6);
  Serial.print(",");
  Serial.print(BENCH_AVG_POWER_MW, 3);
  Serial.print(",");
  Serial.print(totalEnergyMJ, 6);
  Serial.print(",");
  Serial.print(energyPerHashMJ, 9);
  Serial.print(",");
  Serial.print(energyPerByteUJ, 9);
  Serial.print(",");
  Serial.print(ESP.getSketchSize());
  Serial.print(",");
  Serial.print(heapBefore);
  Serial.print(",");
  Serial.print(heapAfter);
  Serial.print(",");
  Serial.print((int32_t)heapBefore - (int32_t)heapAfter);
  Serial.print(",");
  Serial.print(minHeapBefore);
  Serial.print(",");
  Serial.print(minHeapAfter);
  Serial.print(",");
  Serial.print(maxAllocBefore);
  Serial.print(",");
  Serial.print(maxAllocAfter);
  Serial.print(",");
  printDigestPrefixHex(bench_digest, 8);
  Serial.println();
}

static void runSizeBenchmark() {
  Serial.println("algorithm,board,cpu_mhz,msg_size_bytes,iterations,elapsed_seconds,avg_time_us,cycles_per_hash,cycles_per_byte,throughput_kb_s,avg_power_mw,total_energy_mJ,energy_per_hash_mJ,energy_per_byte_uJ,rom_sketch_bytes,heap_before_bytes,heap_after_bytes,heap_delta_bytes,min_heap_before_bytes,min_heap_after_bytes,max_alloc_before_bytes,max_alloc_after_bytes,digest_prefix64_hex");

  for (size_t i = 0; i < sizeof(BENCH_SIZES) / sizeof(BENCH_SIZES[0]); ++i) {
    runOneSizeBenchmark(BENCH_SIZES[i]);
    delay(200);
  }
}

static void runEnduranceTest() {
  fillPayload(bench_payload, ENDURANCE_PAYLOAD_BYTES, 0xC001D00DUL);

  Serial.println();
  Serial.println("ENDURANCE_STATUS,ARMING_10_SECONDS");
  Serial.println("Reset or start FNIRSI FNB58 logging now. The measured hashing phase starts after this countdown.");
  Serial.flush();

  for (int countdown = 10; countdown > 0; --countdown) {
    Serial.print("ENDURANCE_STARTS_IN_SECONDS,");
    Serial.println(countdown);
    Serial.flush();
    delay(1000);
  }

  Serial.println("ENDURANCE_STATUS,START");
  Serial.println("Use FNIRSI FNB58 to record average power and total energy for this phase.");
  Serial.flush();

  const uint64_t startUs = benchNowUs();
  const uint64_t durationUs = (uint64_t)ENDURANCE_SECONDS * 1000000ULL;
  uint64_t nextProgressUs = (uint64_t)ENDURANCE_PROGRESS_SECONDS * 1000000ULL;
  uint64_t totalHashes = 0;

  while ((benchNowUs() - startUs) < durationUs) {
    bench_payload[0] ^= static_cast<uint8_t>(totalHashes & 0xFFU);
    runHashFunction(bench_payload, ENDURANCE_PAYLOAD_BYTES, bench_digest);
    bench_payload[0] ^= static_cast<uint8_t>(totalHashes & 0xFFU);
    ++totalHashes;

    const uint64_t elapsedUs = benchNowUs() - startUs;
    if (elapsedUs >= nextProgressUs) {
      Serial.print("ENDURANCE_PROGRESS,seconds=");
      Serial.print(elapsedUs / 1000000ULL);
      Serial.print(",hashes=");
      Serial.println((unsigned long long)totalHashes);
      Serial.flush();
      nextProgressUs += (uint64_t)ENDURANCE_PROGRESS_SECONDS * 1000000ULL;
    }
  }

  const uint64_t endUs = benchNowUs();
  const double elapsedSec = static_cast<double>(endUs - startUs) / 1000000.0;
  const double avgTimeUs = (elapsedSec * 1000000.0) / static_cast<double>(totalHashes);
  const double totalBytes = static_cast<double>(totalHashes) * static_cast<double>(ENDURANCE_PAYLOAD_BYTES);
  const double throughputKBs = totalBytes / elapsedSec / 1024.0;
  const double cyclesPerHash = avgTimeUs * static_cast<double>(ESP.getCpuFreqMHz());
  const double cyclesPerByte = cyclesPerHash / static_cast<double>(ENDURANCE_PAYLOAD_BYTES);
  const double totalEnergyMJ = (BENCH_AVG_POWER_MW > 0.0) ? (BENCH_AVG_POWER_MW * elapsedSec) : 0.0;
  const double energyPerHashMJ = (BENCH_AVG_POWER_MW > 0.0) ? (totalEnergyMJ / static_cast<double>(totalHashes)) : 0.0;
  const double energyPerByteUJ = (BENCH_AVG_POWER_MW > 0.0) ? ((totalEnergyMJ * 1000.0) / totalBytes) : 0.0;

  Serial.println();
  Serial.println("algorithm,board,cpu_mhz,payload_bytes,duration_seconds,total_hashes,total_bytes,avg_time_us,cycles_per_hash,cycles_per_byte,throughput_kb_s,avg_power_mw,total_energy_mJ,energy_per_hash_mJ,energy_per_byte_uJ,rom_sketch_bytes,free_heap_bytes,min_free_heap_bytes,max_alloc_heap_bytes,digest_prefix64_hex");
  Serial.print(BENCH_ALGORITHM);
  Serial.print(",ESP32,");
  Serial.print(ESP.getCpuFreqMHz());
  Serial.print(",");
  Serial.print(ENDURANCE_PAYLOAD_BYTES);
  Serial.print(",");
  Serial.print(elapsedSec, 6);
  Serial.print(",");
  Serial.print((unsigned long long)totalHashes);
  Serial.print(",");
  Serial.print((unsigned long long)(totalHashes * ENDURANCE_PAYLOAD_BYTES));
  Serial.print(",");
  Serial.print(avgTimeUs, 3);
  Serial.print(",");
  Serial.print(cyclesPerHash, 3);
  Serial.print(",");
  Serial.print(cyclesPerByte, 3);
  Serial.print(",");
  Serial.print(throughputKBs, 6);
  Serial.print(",");
  Serial.print(BENCH_AVG_POWER_MW, 3);
  Serial.print(",");
  Serial.print(totalEnergyMJ, 6);
  Serial.print(",");
  Serial.print(energyPerHashMJ, 9);
  Serial.print(",");
  Serial.print(energyPerByteUJ, 9);
  Serial.print(",");
  Serial.print(ESP.getSketchSize());
  Serial.print(",");
  Serial.print(ESP.getFreeHeap());
  Serial.print(",");
  Serial.print(ESP.getMinFreeHeap());
  Serial.print(",");
  Serial.print(ESP.getMaxAllocHeap());
  Serial.print(",");
  printDigestPrefixHex(bench_digest, 8);
  Serial.println();
  Serial.println("ENDURANCE_STATUS,END");
}

void setup() {
  Serial.begin(BENCH_SERIAL_BAUD);
  delay(1500);

  setCpuFrequencyMhz(BENCH_CPU_MHZ);
  disableRadiosForMeasurement();
  delay(500);

  printSystemMetadata();

  runHashFunction((const uint8_t*)"abc", 3, bench_digest);
  Serial.print("kat_abc_digest_prefix64_hex=");
  printDigestPrefixHex(bench_digest, 8);
  Serial.println();

  if (RUN_SIZE_BENCHMARK) {
    runSizeBenchmark();
  }

  if (RUN_ENDURANCE_TEST) {
    runEnduranceTest();
  } else {
    Serial.println();
    Serial.println("ENDURANCE_STATUS,DISABLED");
  }

  Serial.println("BENCHMARK_STATUS,DONE");
}

void loop() {
  delay(1000);
}
