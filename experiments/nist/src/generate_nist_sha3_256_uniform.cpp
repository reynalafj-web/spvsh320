#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

// ============================================================
//  UNIFORM NIST STS DATA GENERATOR: SHA-3-256
// ============================================================
//
//  Methodology:
//  - Candidate input      : 8-byte little-endian counter LE64(i)
//  - Hash output          : 256-bit digest
//  - Output stream format : ASCII bitstream, characters '0' and '1'
//  - Target size          : 1,000,000,000 bits
//  - NIST STS setting     : 1000 sequences x 1,000,000 bits
//
//  This generator uses a standalone SHA3-256 implementation with
//  Keccak-f[1600], 24 rounds, SHA-3 domain suffix 0x06, and 256-bit output.
// ============================================================

static inline uint64_t rotl64(uint64_t x, unsigned y) {
    return (x << y) | (x >> (64 - y));
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

static const int keccakf_rotc[24] = {
    1, 3, 6, 10, 15, 21, 28, 36, 45, 55, 2, 14,
    27, 41, 56, 8, 25, 43, 62, 18, 39, 61, 20, 44
};

static const int keccakf_piln[24] = {
    10, 7, 11, 17, 18, 3, 5, 16, 8, 21, 24, 4,
    15, 23, 19, 13, 12, 2, 20, 14, 22, 9, 6, 1
};

static void keccakf(uint64_t st[25]) {
    uint64_t bc[5];

    for (int round = 0; round < 24; ++round) {
        for (int i = 0; i < 5; ++i) {
            bc[i] = st[i] ^ st[i + 5] ^ st[i + 10] ^ st[i + 15] ^ st[i + 20];
        }

        for (int i = 0; i < 5; ++i) {
            uint64_t t = bc[(i + 4) % 5] ^ rotl64(bc[(i + 1) % 5], 1);

            for (int j = 0; j < 25; j += 5) {
                st[j + i] ^= t;
            }
        }

        uint64_t t = st[1];

        for (int i = 0; i < 24; ++i) {
            int j = keccakf_piln[i];
            bc[0] = st[j];
            st[j] = rotl64(t, keccakf_rotc[i]);
            t = bc[0];
        }

        for (int j = 0; j < 25; j += 5) {
            for (int i = 0; i < 5; ++i) {
                bc[i] = st[j + i];
            }

            for (int i = 0; i < 5; ++i) {
                st[j + i] ^= (~bc[(i + 1) % 5]) & bc[(i + 2) % 5];
            }
        }

        st[0] ^= keccakf_rndc[round];
    }
}

struct Sha3Ctx {
    union {
        uint8_t b[200];
        uint64_t q[25];
    } st;

    std::size_t pt = 0;
    std::size_t rsiz = 0;
    std::size_t mdlen = 0;
};

static void sha3_256_init(Sha3Ctx& c) {
    std::memset(&c, 0, sizeof(c));
    c.mdlen = 32;
    c.rsiz = 200 - 2 * c.mdlen;
    c.pt = 0;
}

static void sha3_256_update(Sha3Ctx& c, const uint8_t* data, std::size_t len) {
    for (std::size_t i = 0; i < len; ++i) {
        c.st.b[c.pt++] ^= data[i];

        if (c.pt >= c.rsiz) {
            keccakf(c.st.q);
            c.pt = 0;
        }
    }
}

static void sha3_256_final(Sha3Ctx& c, uint8_t out[32]) {
    c.st.b[c.pt] ^= 0x06;
    c.st.b[c.rsiz - 1] ^= 0x80;
    keccakf(c.st.q);
    std::memcpy(out, c.st.b, 32);
}

static void hash_sha3_256(const uint8_t* input, std::size_t len, uint8_t out[32]) {
    Sha3Ctx ctx;
    sha3_256_init(ctx);
    sha3_256_update(ctx, input, len);
    sha3_256_final(ctx, out);
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

    std::string filename = "sha3_256_nist_1Gbit_ascii.txt";

    if (argc >= 2) {
        filename = argv[1];
    }

    std::cout << "============================================================\n";
    std::cout << "  UNIFORM NIST STS DATA GENERATOR\n";
    std::cout << "  Algorithm  : SHA-3-256\n";
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
        hash_sha3_256(payload, 8, digest);
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
