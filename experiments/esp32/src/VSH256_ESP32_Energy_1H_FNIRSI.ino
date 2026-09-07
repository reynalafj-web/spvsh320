/*
  VSH-256 ESP32 Benchmark for Arduino IDE

  Purpose:
    Measure implementation metrics of the VSH-256 baseline on ESP32.

  Metrics printed through Serial Monitor:
    - average execution time per hash
    - derived cycles per hash at 240 MHz
    - cycles per byte
    - throughput in KB/s
    - heap before/after
    - digest prefix sanity check
    - optional 1-hour endurance counter for FNIRSI FNB58 energy measurement

  Notes:
    - This implementation keeps the VSH-256 construction used in the thesis
      experiments: fixed 256-bit RSA-style modulus n = p*q, small-prime VSH
      compression, zero block padding, binary length appending, and fixed
      256-bit big-endian output serialization.
    - No Boost or external big-integer library is used. A minimal 256-bit
      integer implementation is included so the sketch can compile in Arduino IDE.
    - The modular multiplication is intentionally simple and portable
      double-and-add modulo n. It is expected to be slow on ESP32, which is part
      of the baseline cost measurement.
*/
#pragma GCC optimize ("O3,unroll-loops")
#include <Arduino.h>
#include <WiFi.h>

#if defined(ESP32)
#include "esp_timer.h"
#endif

// -----------------------------------------------------------------------------
// Standardized benchmark constants are defined in the benchmark harness below.
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// VSH-256 implementation namespace
// -----------------------------------------------------------------------------
namespace VSH256Impl {

static const size_t LIMBS = 8; // 8 x 32-bit = 256-bit

struct U256 {
  uint32_t w[LIMBS]; // little-endian 32-bit limbs
};

static inline U256 zero256() {
  U256 x;
  for (size_t i = 0; i < LIMBS; ++i) x.w[i] = 0;
  return x;
}

static inline U256 one256() {
  U256 x = zero256();
  x.w[0] = 1;
  return x;
}

static inline bool isZero(const U256& x) {
  for (size_t i = 0; i < LIMBS; ++i) if (x.w[i] != 0) return false;
  return true;
}

static int cmp256(const U256& a, const U256& b) {
  for (int i = static_cast<int>(LIMBS) - 1; i >= 0; --i) {
    if (a.w[i] < b.w[i]) return -1;
    if (a.w[i] > b.w[i]) return 1;
  }
  return 0;
}

static void subInPlace(U256& a, const U256& b) {
  uint64_t borrow = 0;
  for (size_t i = 0; i < LIMBS; ++i) {
    uint64_t ai = static_cast<uint64_t>(a.w[i]);
    uint64_t bi = static_cast<uint64_t>(b.w[i]) + borrow;
    if (ai >= bi) {
      a.w[i] = static_cast<uint32_t>(ai - bi);
      borrow = 0;
    } else {
      a.w[i] = static_cast<uint32_t>((uint64_t(1) << 32) + ai - bi);
      borrow = 1;
    }
  }
}

static uint32_t addRaw(U256& out, const U256& a, const U256& b) {
  uint64_t carry = 0;
  for (size_t i = 0; i < LIMBS; ++i) {
    uint64_t s = static_cast<uint64_t>(a.w[i]) + b.w[i] + carry;
    out.w[i] = static_cast<uint32_t>(s & 0xFFFFFFFFULL);
    carry = s >> 32;
  }
  return static_cast<uint32_t>(carry);
}

// n = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFE4E0000000000000000000000000000AACD
static const U256 MOD_N = {{
  0x0000AACDUL,
  0x00000000UL,
  0x00000000UL,
  0x00000000UL,
  0xFFFFFE4EUL,
  0xFFFFFFFFUL,
  0xFFFFFFFFUL,
  0xFFFFFFFFUL
}};

static U256 addMod(const U256& a, const U256& b) {
  U256 s;
  uint32_t carry = addRaw(s, a, b);

  if (carry || cmp256(s, MOD_N) >= 0) {
    subInPlace(s, MOD_N);
  }
  return s;
}

static U256 doubleMod(const U256& a) {
  return addMod(a, a);
}

static bool getBit(const U256& x, size_t bit_index) {
  size_t limb = bit_index / 32;
  size_t bit = bit_index % 32;
  if (limb >= LIMBS) return false;
  return ((x.w[limb] >> bit) & 1U) != 0;
}

static U256 mulMod(U256 a, const U256& b) {
  U256 result = zero256();

  for (size_t bit = 0; bit < 256; ++bit) {
    if (getBit(b, bit)) {
      result = addMod(result, a);
    }
    a = doubleMod(a);
  }
  return result;
}

static U256 mulSmall(const U256& a, uint32_t m) {
  U256 out;
  uint64_t carry = 0;
  for (size_t i = 0; i < LIMBS; ++i) {
    uint64_t p = static_cast<uint64_t>(a.w[i]) * static_cast<uint64_t>(m) + carry;
    out.w[i] = static_cast<uint32_t>(p & 0xFFFFFFFFULL);
    carry = p >> 32;
  }
  return out;
}

static const uint32_t VSH_PRIMES[43] = {
  2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
  37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
  79, 83, 89, 97, 101, 103, 107, 109, 113,
  127, 131, 137, 139, 149, 151, 157, 163,
  167, 173, 179, 181, 191
};

static const size_t VSH_K = 43;

static bool getMessageBitMSB(const uint8_t* msg, size_t len_bytes, size_t bit_index) {
  size_t byte_index = bit_index / 8;
  size_t bit_in_byte = bit_index % 8;
  if (byte_index >= len_bytes) return false;
  return ((msg[byte_index] >> (7 - bit_in_byte)) & 1U) != 0;
}

static size_t bitLengthSizeT(size_t x) {
  if (x == 0) return 1;
  size_t bits = 0;
  while (x != 0) {
    ++bits;
    x >>= 1;
  }
  return bits;
}

static bool getLengthBitMSB(size_t value, size_t bit_index, size_t total_bits) {
  if (total_bits == 0) return false;
  size_t shift = total_bits - 1 - bit_index;
  return ((value >> shift) & 1U) != 0;
}

static bool getVSHInputBit(const uint8_t* msg,
                           size_t len_bytes,
                           size_t bit_index,
                           size_t msg_bits,
                           size_t padded_len,
                           size_t len_bits) {
  if (bit_index < msg_bits) {
    return getMessageBitMSB(msg, len_bytes, bit_index);
  }
  if (bit_index < padded_len) {
    return false; // zero padding to k-bit block boundary
  }
  if (bit_index < padded_len + len_bits) {
    return getLengthBitMSB(msg_bits, bit_index - padded_len, len_bits);
  }
  return false;
}

static U256 VSHHash256(const uint8_t* msg, size_t len_bytes) {
  const size_t msg_bits = len_bytes * 8;
  const size_t L = (msg_bits + VSH_K - 1) / VSH_K;
  const size_t padded_len = L * VSH_K;
  const size_t len_bits = bitLengthSizeT(msg_bits);
  const size_t total_bits = padded_len + len_bits;
  const size_t total_blocks = (total_bits + VSH_K - 1) / VSH_K;

  U256 x = one256();

  for (size_t block = 0; block < total_blocks; ++block) {
    U256 acc = one256();

    for (size_t i = 0; i < VSH_K; ++i) {
      size_t bit_index = block * VSH_K + i;
      bool bit = getVSHInputBit(msg, len_bytes, bit_index, msg_bits, padded_len, len_bits);
      if (bit) {
        acc = mulSmall(acc, VSH_PRIMES[i]);
      }
    }

    x = mulMod(x, x);
    x = mulMod(x, acc);
  }

  return x;
}

static void U256ToBytesBE(const U256& x, uint8_t out[32]) {
  // U256 limbs are little-endian 32-bit words. Output is fixed 256-bit big-endian.
  for (int limb = 7; limb >= 0; --limb) {
    uint32_t w = x.w[limb];
    size_t offset = static_cast<size_t>(7 - limb) * 4;
    out[offset + 0] = static_cast<uint8_t>((w >> 24) & 0xFFU);
    out[offset + 1] = static_cast<uint8_t>((w >> 16) & 0xFFU);
    out[offset + 2] = static_cast<uint8_t>((w >> 8) & 0xFFU);
    out[offset + 3] = static_cast<uint8_t>(w & 0xFFU);
  }
}

static void VSH256_Hash(const uint8_t* msg, size_t len_bytes, uint8_t out[32]) {
  U256 h = VSHHash256(msg, len_bytes);
  U256ToBytesBE(h, out);
}

} // namespace VSH256Impl


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

static const char* BENCH_ALGORITHM = "VSH-256";
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
  VSH256Impl::VSH256_Hash(msg, len, out);
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
