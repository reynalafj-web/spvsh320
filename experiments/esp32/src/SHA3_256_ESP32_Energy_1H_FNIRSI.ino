/*
 * SHA3-256 ESP32 Benchmark for Arduino IDE
 *
 * Purpose:
 *   Microcontroller-board implementation metrics for SHA3-256 on ESP32.
 *   The benchmark reports execution time, derived cycles, cycles per byte,
 *   throughput, and heap readings. It also runs a 1-hour
 *   endurance test for FNIRSI FNB58 board-level power/energy measurement.
 *
 * Notes:
 *   - This file is written as a single .ino sketch.
 *   - The SHA3 context type is declared before all SHA3 functions to avoid
 *     Arduino IDE auto-prototype issues.
 *   - SHA3-256 uses the standard SHA-3 domain suffix 0x06.
 *   - Wi-Fi is disabled in setup for more stable board-level measurement.
 */
#pragma GCC optimize ("O3,unroll-loops")
#include <Arduino.h>
#include "esp_timer.h"
#include <WiFi.h>
#include <cstring>
#include <cstdint>

// ============================================================
// SHA3-256 core implementation
// ============================================================
static const size_t SHA3_STATE_WORDS = 25;
static const size_t SHA3_256_RATE_BYTES = 136; // 1088-bit rate
static const int SHA3_KECCAK_ROUNDS = 24;

typedef struct {
  uint64_t s[SHA3_STATE_WORDS];
} sha3_context;

static inline uint64_t rotl64(uint64_t x, unsigned r) {
  return (x << r) | (x >> (64 - r));
}

static inline uint64_t load64_le(const uint8_t* p) {
  uint64_t v = 0;
  for (int i = 7; i >= 0; --i) {
    v = (v << 8) | p[i];
  }
  return v;
}

static inline void store64_le(uint64_t v, uint8_t* p) {
  for (int i = 0; i < 8; ++i) {
    p[i] = static_cast<uint8_t>(v & 0xFFU);
    v >>= 8;
  }
}

static const uint64_t keccakf_rndc[24] = {
  UINT64_C(0x0000000000000001), UINT64_C(0x0000000000008082),
  UINT64_C(0x800000000000808a), UINT64_C(0x8000000080008000),
  UINT64_C(0x000000000000808b), UINT64_C(0x0000000080000001),
  UINT64_C(0x8000000080008081), UINT64_C(0x8000000000008009),
  UINT64_C(0x000000000000008a), UINT64_C(0x0000000000000088),
  UINT64_C(0x0000000080008009), UINT64_C(0x000000008000000a),
  UINT64_C(0x000000008000808b), UINT64_C(0x800000000000008b),
  UINT64_C(0x8000000000008089), UINT64_C(0x8000000000008003),
  UINT64_C(0x8000000000008002), UINT64_C(0x8000000000000080),
  UINT64_C(0x000000000000800a), UINT64_C(0x800000008000000a),
  UINT64_C(0x8000000080008081), UINT64_C(0x8000000000008080),
  UINT64_C(0x0000000080000001), UINT64_C(0x8000000080008008)
};

static const unsigned keccakf_rotc[24] = {
  1, 3, 6, 10, 15, 21, 28, 36, 45, 55, 2, 14,
  27, 41, 56, 8, 25, 43, 62, 18, 39, 61, 20, 44
};

static const unsigned keccakf_piln[24] = {
  10, 7, 11, 17, 18, 3, 5, 16, 8, 21, 24, 4,
  15, 23, 19, 13, 12, 2, 20, 14, 22, 9, 6, 1
};

static void keccakf(uint64_t st[25]) {
  uint64_t bc[5];

  for (int round = 0; round < SHA3_KECCAK_ROUNDS; ++round) {
    // Theta
    for (int i = 0; i < 5; ++i) {
      bc[i] = st[i] ^ st[i + 5] ^ st[i + 10] ^ st[i + 15] ^ st[i + 20];
    }

    for (int i = 0; i < 5; ++i) {
      uint64_t t = bc[(i + 4) % 5] ^ rotl64(bc[(i + 1) % 5], 1);
      for (int j = 0; j < 25; j += 5) {
        st[j + i] ^= t;
      }
    }

    // Rho and Pi
    uint64_t t = st[1];
    for (int i = 0; i < 24; ++i) {
      int j = keccakf_piln[i];
      uint64_t tmp = st[j];
      st[j] = rotl64(t, keccakf_rotc[i]);
      t = tmp;
    }

    // Chi
    for (int j = 0; j < 25; j += 5) {
      for (int i = 0; i < 5; ++i) {
        bc[i] = st[j + i];
      }
      for (int i = 0; i < 5; ++i) {
        st[j + i] ^= (~bc[(i + 1) % 5]) & bc[(i + 2) % 5];
      }
    }

    // Iota
    st[0] ^= keccakf_rndc[round];
  }
}

static void sha3_256_hash(const uint8_t* in, size_t inlen, uint8_t out[32]) {
  sha3_context ctx;
  std::memset(&ctx, 0, sizeof(ctx));

  // Absorb full rate blocks.
  while (inlen >= SHA3_256_RATE_BYTES) {
    for (size_t i = 0; i < SHA3_256_RATE_BYTES / 8; ++i) {
      ctx.s[i] ^= load64_le(in + (8 * i));
    }
    keccakf(ctx.s);
    in += SHA3_256_RATE_BYTES;
    inlen -= SHA3_256_RATE_BYTES;
  }

  // Final block with SHA-3 domain suffix 0x06 and pad10*1.
  uint8_t block[SHA3_256_RATE_BYTES];
  std::memset(block, 0, sizeof(block));
  if (inlen > 0) {
    std::memcpy(block, in, inlen);
  }
  block[inlen] ^= 0x06;
  block[SHA3_256_RATE_BYTES - 1] ^= 0x80;

  for (size_t i = 0; i < SHA3_256_RATE_BYTES / 8; ++i) {
    ctx.s[i] ^= load64_le(block + (8 * i));
  }
  keccakf(ctx.s);

  // Squeeze 32 bytes.
  for (size_t i = 0; i < 4; ++i) {
    store64_le(ctx.s[i], out + (8 * i));
  }
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

static const char* BENCH_ALGORITHM = "SHA3-256";
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
  sha3_256_hash(msg, len, out);
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
