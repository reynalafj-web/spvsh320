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
// Benchmark configuration
// -----------------------------------------------------------------------------
static const uint32_t CPU_CLOCK_MHZ = 240;
static const bool RUN_ENDURANCE_TEST = false;   // set true for 1-hour FNB58 run
static const uint32_t ENDURANCE_SECONDS = 3600;
static const size_t ENDURANCE_PAYLOAD_BYTES = 1024;
static const uint32_t SERIAL_BAUD = 115200;

struct TestCase {
  size_t message_size;
  uint32_t iterations;
};

// VSH-256 is intentionally much slower than SpVSH-320 and SHA3-256 on ESP32.
// Keep default iterations small, then increase manually if the run is too short.
static const TestCase TEST_CASES[] = {
  {16,   10},
  {128,   5},
  {242,   3},
  {1024,  1}
};

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

// -----------------------------------------------------------------------------
// Benchmark utilities
// -----------------------------------------------------------------------------
static void fillPayload(uint8_t* buf, size_t len, uint32_t seed) {
  uint32_t x = seed;
  for (size_t i = 0; i < len; ++i) {
    x = x * 1664525UL + 1013904223UL;
    buf[i] = static_cast<uint8_t>((x >> 24) & 0xFFU);
  }
}

static void printDigestPrefix(const uint8_t digest[32]) {
  for (int i = 0; i < 8; ++i) {
    if (digest[i] < 16) Serial.print('0');
    Serial.print(digest[i], HEX);
  }
}

static uint64_t nowUs() {
#if defined(ESP32)
  return static_cast<uint64_t>(esp_timer_get_time());
#else
  return static_cast<uint64_t>(micros());
#endif
}

static void runLatencyBenchmark() {
  static uint8_t payload[1024];
  uint8_t digest[32];

  Serial.println();
  Serial.println("============================================================");
  Serial.println("  VSH-256 ESP32 IMPLEMENTATION BENCHMARK");
  Serial.println("============================================================");
  Serial.println("algorithm,board,cpu_clock_mhz,message_size_bytes,iterations,total_time_us,avg_time_us,cycles_per_hash,cycles_per_byte,throughput_kb_s,heap_before,heap_after,heap_used,digest_prefix");

  const size_t test_count = sizeof(TEST_CASES) / sizeof(TEST_CASES[0]);

  for (size_t tc = 0; tc < test_count; ++tc) {
    size_t msg_len = TEST_CASES[tc].message_size;
    uint32_t iterations = TEST_CASES[tc].iterations;

    fillPayload(payload, msg_len, 0x2560BEEFUL + static_cast<uint32_t>(msg_len));

    // Warm-up hash outside measured loop.
    VSH256Impl::VSH256_Hash(payload, msg_len, digest);

    uint32_t heap_before = ESP.getFreeHeap();
    uint64_t start = nowUs();

    for (uint32_t i = 0; i < iterations; ++i) {
      payload[0] ^= static_cast<uint8_t>(i & 0xFFU); // prevent over-aggressive optimization
      VSH256Impl::VSH256_Hash(payload, msg_len, digest);
      payload[0] ^= static_cast<uint8_t>(i & 0xFFU);
    }

    uint64_t end = nowUs();
    uint32_t heap_after = ESP.getFreeHeap();

    double total_us = static_cast<double>(end - start);
    double avg_us = total_us / static_cast<double>(iterations);
    double cycles_per_hash = avg_us * static_cast<double>(CPU_CLOCK_MHZ);
    double cycles_per_byte = cycles_per_hash / static_cast<double>(msg_len);
    double throughput_kb_s = (static_cast<double>(msg_len) * static_cast<double>(iterations) * 1000000.0) / (total_us * 1024.0);
    int32_t heap_used = static_cast<int32_t>(heap_before) - static_cast<int32_t>(heap_after);

    Serial.print("VSH-256,ESP32 DevKit,");
    Serial.print(CPU_CLOCK_MHZ);
    Serial.print(',');
    Serial.print(msg_len);
    Serial.print(',');
    Serial.print(iterations);
    Serial.print(',');
    Serial.print(total_us, 2);
    Serial.print(',');
    Serial.print(avg_us, 2);
    Serial.print(',');
    Serial.print(cycles_per_hash, 2);
    Serial.print(',');
    Serial.print(cycles_per_byte, 2);
    Serial.print(',');
    Serial.print(throughput_kb_s, 4);
    Serial.print(',');
    Serial.print(heap_before);
    Serial.print(',');
    Serial.print(heap_after);
    Serial.print(',');
    Serial.print(heap_used);
    Serial.print(',');
    printDigestPrefix(digest);
    Serial.println();
  }
}

static void runEnduranceTest() {
  static uint8_t payload1024[ENDURANCE_PAYLOAD_BYTES];
  uint8_t digest[32];

  fillPayload(payload1024, ENDURANCE_PAYLOAD_BYTES, 0x2560B58UL);

  Serial.println();
  Serial.println("============================================================");
  Serial.println("  VSH-256 1-HOUR ENDURANCE TEST");
  Serial.println("============================================================");
  Serial.println("Use FNIRSI FNB58 to record active energy during this phase.");
  Serial.println("Avoid Serial Monitor polling if you want the cleanest power reading.");
  Serial.println();
  Serial.print("payload_size_bytes=");
  Serial.println(ENDURANCE_PAYLOAD_BYTES);
  Serial.print("target_duration_seconds=");
  Serial.println(ENDURANCE_SECONDS);
  Serial.println("status=START");
  Serial.flush();

  uint64_t start_us = nowUs();
  uint64_t duration_us = static_cast<uint64_t>(ENDURANCE_SECONDS) * 1000000ULL;
  uint64_t total_hashes = 0;

  while ((nowUs() - start_us) < duration_us) {
    payload1024[0] ^= static_cast<uint8_t>(total_hashes & 0xFFU);
    VSH256Impl::VSH256_Hash(payload1024, ENDURANCE_PAYLOAD_BYTES, digest);
    payload1024[0] ^= static_cast<uint8_t>(total_hashes & 0xFFU);
    ++total_hashes;
  }

  uint64_t end_us = nowUs();
  double elapsed_seconds = static_cast<double>(end_us - start_us) / 1000000.0;
  double avg_time_us = (elapsed_seconds * 1000000.0) / static_cast<double>(total_hashes);
  double throughput_kb_s = (static_cast<double>(total_hashes) * static_cast<double>(ENDURANCE_PAYLOAD_BYTES)) / (elapsed_seconds * 1024.0);
  double cycles_per_hash = avg_time_us * static_cast<double>(CPU_CLOCK_MHZ);
  double cycles_per_byte = cycles_per_hash / static_cast<double>(ENDURANCE_PAYLOAD_BYTES);

  Serial.println("status=END");
  Serial.print("algorithm=VSH-256\n");
  Serial.print("duration_seconds=");
  Serial.println(elapsed_seconds, 3);
  Serial.print("total_hashes_completed=");
  Serial.println(static_cast<unsigned long long>(total_hashes));
  Serial.print("total_bytes_processed=");
  Serial.println(static_cast<unsigned long long>(total_hashes * ENDURANCE_PAYLOAD_BYTES));
  Serial.print("avg_time_us=");
  Serial.println(avg_time_us, 3);
  Serial.print("cycles_per_hash=");
  Serial.println(cycles_per_hash, 3);
  Serial.print("cycles_per_byte=");
  Serial.println(cycles_per_byte, 3);
  Serial.print("throughput_kb_s=");
  Serial.println(throughput_kb_s, 6);
  Serial.print("digest_prefix=");
  printDigestPrefix(digest);
  Serial.println();
}

void setup() {
  Serial.begin(SERIAL_BAUD);
  delay(1500);

#if defined(ESP32)
  setCpuFrequencyMhz(CPU_CLOCK_MHZ);
#endif

  WiFi.mode(WIFI_OFF);
  btStop();

  Serial.println();
  Serial.println("VSH-256 ESP32 Benchmark");
  Serial.println("Wi-Fi and Bluetooth disabled.");
  Serial.print("CPU frequency MHz: ");
  Serial.println(getCpuFrequencyMhz());
  Serial.print("Free heap at start: ");
  Serial.println(ESP.getFreeHeap());

  runLatencyBenchmark();

  if (RUN_ENDURANCE_TEST) {
    runEnduranceTest();
  } else {
    Serial.println();
    Serial.println("Endurance test is disabled. Set RUN_ENDURANCE_TEST=true for the 1-hour FNB58 run.");
  }
}

void loop() {
  delay(1000);
}
