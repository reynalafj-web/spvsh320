/*
 * SHA3-256 ESP32 Benchmark for Arduino IDE
 *
 * Purpose:
 *   Microcontroller-board implementation metrics for SHA3-256 on ESP32.
 *   The benchmark reports execution time, derived cycles, cycles per byte,
 *   throughput, and heap readings. It also supports an optional 1-hour
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
#include <WiFi.h>
#include <cstring>
#include <cstdint>

// ============================================================
// Benchmark configuration
// ============================================================
static const uint32_t CPU_MHZ = 240;
static const uint32_t CPU_HZ  = CPU_MHZ * 1000000UL;

static const size_t DIGEST_BYTES = 32;
static const size_t MAX_PAYLOAD_BYTES = 1024;

// Change this to true for the 1-hour FNIRSI FNB58 endurance test.
static const bool RUN_ENDURANCE_TEST = false;
static const uint32_t ENDURANCE_SECONDS = 3600;
static const size_t ENDURANCE_PAYLOAD_BYTES = 1024;

// Iterations for short benchmark runs.
// Increase if your timing is noisy.
static const uint32_t DEFAULT_ITERATIONS = 1000;

static uint8_t payload[MAX_PAYLOAD_BYTES];
static uint8_t digest[DIGEST_BYTES];

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
// Benchmark helpers
// ============================================================
static void fillPayload(uint8_t* buf, size_t len, uint32_t seed) {
  uint32_t x = seed;
  for (size_t i = 0; i < len; ++i) {
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    buf[i] = static_cast<uint8_t>(x & 0xFFU);
  }
}

static String digestPrefixHex(const uint8_t* d, size_t prefixBytes) {
  static const char* hex = "0123456789abcdef";
  String s;
  for (size_t i = 0; i < prefixBytes; ++i) {
    s += hex[(d[i] >> 4) & 0xF];
    s += hex[d[i] & 0xF];
  }
  return s;
}

static void runLatencyBenchmark(size_t msgLen, uint32_t iterations) {
  fillPayload(payload, msgLen, 0x5A3A0256UL + static_cast<uint32_t>(msgLen));

  // Warm-up phase.
  for (uint32_t i = 0; i < 50; ++i) {
    sha3_256_hash(payload, msgLen, digest);
  }

  uint32_t heapBefore = ESP.getFreeHeap();
  int64_t startUs = esp_timer_get_time();

  for (uint32_t i = 0; i < iterations; ++i) {
    // Change one byte very lightly so the compiler cannot treat the input as constant.
    payload[0] = static_cast<uint8_t>(payload[0] + 1);
    sha3_256_hash(payload, msgLen, digest);
  }

  int64_t endUs = esp_timer_get_time();
  uint32_t heapAfter = ESP.getFreeHeap();

  double totalUs = static_cast<double>(endUs - startUs);
  double avgUs = totalUs / static_cast<double>(iterations);
  double cyclesPerHash = avgUs * static_cast<double>(CPU_MHZ);
  double cyclesPerByte = cyclesPerHash / static_cast<double>(msgLen == 0 ? 1 : msgLen);
  double throughputKBs = (static_cast<double>(msgLen) * static_cast<double>(iterations) * 1000000.0) /
                         (totalUs * 1024.0);

  Serial.print("SHA3-256,ESP32,");
  Serial.print(CPU_MHZ);
  Serial.print(",");
  Serial.print(msgLen);
  Serial.print(",");
  Serial.print(iterations);
  Serial.print(",");
  Serial.print(avgUs, 4);
  Serial.print(",");
  Serial.print(cyclesPerHash, 2);
  Serial.print(",");
  Serial.print(cyclesPerByte, 2);
  Serial.print(",");
  Serial.print(throughputKBs, 4);
  Serial.print(",");
  Serial.print(heapBefore);
  Serial.print(",");
  Serial.print(heapAfter);
  Serial.print(",");
  Serial.println(digestPrefixHex(digest, 8));
}

static void runEnduranceTest() {
  fillPayload(payload, ENDURANCE_PAYLOAD_BYTES, 0x320F0B58UL);

  Serial.println();
  Serial.println("============================================================");
  Serial.println("SHA3-256 ESP32 1-hour endurance test");
  Serial.println("Use FNIRSI FNB58 to record board-level voltage/current/power/energy.");
  Serial.println("Avoid opening Serial Monitor repeatedly during the active measurement.");
  Serial.println("============================================================");
  Serial.flush();

  const int64_t durationUs = static_cast<int64_t>(ENDURANCE_SECONDS) * 1000000LL;
  const int64_t startUs = esp_timer_get_time();
  int64_t nextProgressUs = startUs + 60000000LL; // every 60 seconds
  uint64_t hashes = 0;

  while ((esp_timer_get_time() - startUs) < durationUs) {
    payload[0] = static_cast<uint8_t>(payload[0] + 1);
    sha3_256_hash(payload, ENDURANCE_PAYLOAD_BYTES, digest);
    ++hashes;

    int64_t now = esp_timer_get_time();
    if (now >= nextProgressUs) {
      uint32_t elapsed = static_cast<uint32_t>((now - startUs) / 1000000LL);
      Serial.print("elapsed_seconds=");
      Serial.print(elapsed);
      Serial.print(", total_hashes=");
      Serial.println(static_cast<unsigned long long>(hashes));
      Serial.flush();
      nextProgressUs += 60000000LL;
    }
  }

  int64_t endUs = esp_timer_get_time();
  double elapsedSec = static_cast<double>(endUs - startUs) / 1000000.0;
  double avgTimeUs = (elapsedSec * 1000000.0) / static_cast<double>(hashes);
  double totalBytes = static_cast<double>(hashes) * static_cast<double>(ENDURANCE_PAYLOAD_BYTES);
  double throughputKBs = totalBytes / elapsedSec / 1024.0;
  double cyclesPerHash = avgTimeUs * static_cast<double>(CPU_MHZ);
  double cyclesPerByte = cyclesPerHash / static_cast<double>(ENDURANCE_PAYLOAD_BYTES);

  Serial.println();
  Serial.println("ENDURANCE_RESULT");
  Serial.print("algorithm=SHA3-256");
  Serial.print(",payload_bytes=");
  Serial.print(ENDURANCE_PAYLOAD_BYTES);
  Serial.print(",elapsed_seconds=");
  Serial.print(elapsedSec, 3);
  Serial.print(",total_hashes=");
  Serial.print(static_cast<unsigned long long>(hashes));
  Serial.print(",total_bytes=");
  Serial.print(static_cast<unsigned long long>(hashes * ENDURANCE_PAYLOAD_BYTES));
  Serial.print(",avg_time_us=");
  Serial.print(avgTimeUs, 4);
  Serial.print(",cycles_per_hash=");
  Serial.print(cyclesPerHash, 2);
  Serial.print(",cycles_per_byte=");
  Serial.print(cyclesPerByte, 2);
  Serial.print(",throughput_kb_s=");
  Serial.print(throughputKBs, 4);
  Serial.print(",digest_prefix=");
  Serial.println(digestPrefixHex(digest, 8));
  Serial.flush();
}

// ============================================================
// Arduino entry points
// ============================================================
void setup() {
  Serial.begin(115200);
  delay(1000);

  setCpuFrequencyMhz(CPU_MHZ);
  WiFi.mode(WIFI_OFF);
  WiFi.disconnect(true);

  Serial.println();
  Serial.println("============================================================");
  Serial.println("SHA3-256 ESP32 Benchmark");
  Serial.println("============================================================");
  Serial.print("CPU MHz: ");
  Serial.println(getCpuFrequencyMhz());
  Serial.print("Free heap: ");
  Serial.println(ESP.getFreeHeap());
  Serial.println();

  Serial.println("algorithm,board,cpu_mhz,msg_size_bytes,iterations,avg_time_us,cycles_per_hash,cycles_per_byte,throughput_kb_s,heap_before,heap_after,digest_prefix64");

  runLatencyBenchmark(16, DEFAULT_ITERATIONS);
  runLatencyBenchmark(64, DEFAULT_ITERATIONS);
  runLatencyBenchmark(128, DEFAULT_ITERATIONS);
  runLatencyBenchmark(242, DEFAULT_ITERATIONS);
  runLatencyBenchmark(256, DEFAULT_ITERATIONS);
  runLatencyBenchmark(512, DEFAULT_ITERATIONS);
  runLatencyBenchmark(1024, DEFAULT_ITERATIONS);

  if (RUN_ENDURANCE_TEST) {
    runEnduranceTest();
  } else {
    Serial.println();
    Serial.println("Endurance test disabled. Set RUN_ENDURANCE_TEST=true for 1-hour FNB58 measurement.");
  }
}

void loop() {
  delay(1000);
}
