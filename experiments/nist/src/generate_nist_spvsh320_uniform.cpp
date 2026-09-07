#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>

// ============================================================
//  UNIFORM NIST STS DATA GENERATOR: SpVSH-320
// ============================================================
//
//  Methodology:
//  - Candidate input      : 8-byte little-endian counter LE64(i)
//  - Hash output          : 256-bit digest
//  - Output stream format : ASCII bitstream, characters '0' and '1'
//  - Target size          : 1,000,000,000 bits
//  - NIST STS setting     : 1000 sequences x 1,000,000 bits
//
//  This generator uses the same SpVSH_StrictSponge_256 function used in
//  the truncated-preimage and truncated-collision scripts.
// ============================================================

struct PentaState {
    uint64_t s[5];
};

const uint64_t P1 = 0xFFFFFFFFFFFFFFC5ULL;
const uint64_t C1 = 59;
const uint64_t P2 = 0xFFFFFFFFFFFFFFADULL;
const uint64_t C2 = 83;

static const uint16_t SMALL_PRIMES[64] = {
      3,   5,   7,  11,  13,  17,  19,  23,  29,  31,  37,  41,  43,  47,  53,  59,
     61,  67,  71,  73,  79,  83,  89,  97, 101, 103, 107, 109, 113, 127, 131, 137,
    139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227,
    229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313
};

static inline void mul64x64_128(uint64_t a, uint64_t b, uint64_t* res_hi, uint64_t* res_lo) {
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

static inline void PermuteSpVSH(PentaState& state) {
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

    state.s[0] = t0 ^ (t1 >> 17) ^ (t1 << 47);
    state.s[1] = t1 ^ (t2 >> 11) ^ (t2 << 53);
    state.s[2] = t2 ^ (t3 >> 19) ^ (t3 << 45);
    state.s[3] = t3 ^ (t4 >> 13) ^ (t4 << 51);
    state.s[4] = t4 ^ (t0 >> 23) ^ (t0 << 41);
}

static inline uint64_t GET_U64(const uint8_t* ptr) {
    uint64_t val;
    std::memcpy(&val, ptr, 8);
    return val;
}

void SpVSH_StrictSponge_256(const uint8_t* in, size_t len, uint8_t* out) {
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

    size_t length = len;

    while (length >= 8) {
        state.s[0] ^= GET_U64(in);
        PermuteSpVSH(state);

        in += 8;
        length -= 8;
    }

    uint8_t buffer[8] = {0};

    if (length > 0) {
        std::memcpy(buffer, in, length);
    }

    buffer[length] = 0x80;
    state.s[0] ^= GET_U64(buffer);
    PermuteSpVSH(state);

    uint64_t h[4];

    for (int i = 0; i < 4; i++) {
        h[i] = state.s[0];

        if (i < 3) {
            PermuteSpVSH(state);
        }
    }

    std::memcpy(out, h, 32);
}

static void counter_to_candidate_le(uint64_t counter, uint8_t payload[8]) {
    for (int i = 0; i < 8; ++i) {
        payload[i] = static_cast<uint8_t>((counter >> (8 * i)) & 0xFFU);
    }
}

static void write_digest_as_ascii_bits(std::ofstream& outfile,
                                       const uint8_t digest[32]) {
    char bits[256];

    for (int b = 0; b < 32; ++b) {
        for (int i = 7; i >= 0; --i) {
            bits[b * 8 + (7 - i)] = ((digest[b] >> i) & 1U) ? '1' : '0';
        }
    }

    outfile.write(bits, 256);
}

int main(int argc, char** argv) {
    const uint64_t TARGET_BITS = 1000000000ULL;
    const int DIGEST_BITS = 256;
    const uint64_t TARGET_ITERATIONS = TARGET_BITS / DIGEST_BITS;

    std::string filename = "spvsh320_nist_1Gbit_ascii.txt";

    if (argc >= 2) {
        filename = argv[1];
    }

    std::cout << "============================================================\n";
    std::cout << "  UNIFORM NIST STS DATA GENERATOR\n";
    std::cout << "  Algorithm  : SpVSH-320\n";
    std::cout << "============================================================\n\n";

    std::cout << "[INFO] Output file        : " << filename << "\n";
    std::cout << "[INFO] Target bits        : " << TARGET_BITS << "\n";
    std::cout << "[INFO] Digest bits        : " << DIGEST_BITS << "\n";
    std::cout << "[INFO] Iterations         : " << TARGET_ITERATIONS << "\n";
    std::cout << "[INFO] Candidate input    : LE64(counter), 8 bytes\n";
    std::cout << "[INFO] Output format      : ASCII bitstream\n";
    std::cout << "[INFO] NIST STS setting   : 1000 sequences x 1,000,000 bits\n\n";

    std::ofstream outfile(filename, std::ios::out | std::ios::binary);

    if (!outfile.is_open()) {
        std::cerr << "[ERROR] Failed to create output file.\n";
        return 1;
    }

    uint8_t payload[8];
    uint8_t digest[32];

    for (uint64_t counter = 0; counter < TARGET_ITERATIONS; ++counter) {
        counter_to_candidate_le(counter, payload);
        SpVSH_StrictSponge_256(payload, 8, digest);
        write_digest_as_ascii_bits(outfile, digest);

        if (counter % 100000 == 0) {
            uint64_t generated_bits = counter * static_cast<uint64_t>(DIGEST_BITS);
            std::cout << "\r[INFO] Progress: " << generated_bits
                      << " / " << TARGET_BITS << " bits" << std::flush;
        }
    }

    outfile.close();

    std::cout << "\r[INFO] Progress: " << TARGET_BITS
              << " / " << TARGET_BITS << " bits\n";
    std::cout << "[DONE] Generated file: " << filename << "\n";

    return 0;
}
