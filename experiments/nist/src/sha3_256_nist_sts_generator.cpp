/*
 * SHA3-256 NIST STS bitstream generator
 *
 * Scenario:
 *   - Generate reproducible bitstreams for NIST SP 800-22 / NIST STS.
 *   - Default dataset: 1000 sequences x 1,000,000 bits = 1,000,000,000 bits.
 *   - Message schedule per digest block:
 *       M_{s,i} = LE64(domain_salt) || LE64(sequence_index) || LE64(block_counter)
 *   - Hash core: SHA3-256, unseeded mode, same SHA3 implementation style used in
 *     the previous SHA-3 experiment script.
 *
 * Output format:
 *   - binary: packed bits, 8 bits per byte, most-significant bit first within each output byte
 *   - ascii : characters '0' and '1'
 */

#include <cassert>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using Clock = std::chrono::high_resolution_clock;
using seed_t = uint64_t;

// ============================================================
// SHA3-256 core implementation
// Kept aligned with the standalone SHA3-256 experiment script.
// ============================================================
#define SHA3_KECCAK_SPONGE_WORDS (((1600) / 8) / sizeof(uint64_t))
#define SHA3_KECCAK_ROUNDS 24

static inline uint64_t ROTL64(uint64_t x, unsigned r) {
    return (x << r) | (x >> (64 - r));
}

template <bool bswap>
static inline uint64_t GET_U64(const uint8_t* p, size_t off) {
    const uint8_t* q = p + off;
    uint64_t v = 0;

    if constexpr (bswap) {
        for (int i = 0; i < 8; ++i) {
            v = (v << 8) | q[i];
        }
    } else {
        for (int i = 7; i >= 0; --i) {
            v = (v << 8) | q[i];
        }
    }

    return v;
}

template <bool bswap>
static inline void PUT_U64(uint64_t v, uint8_t* p, size_t off) {
    uint8_t* q = p + off;

    if constexpr (bswap) {
        for (int i = 7; i >= 0; --i) {
            q[i] = static_cast<uint8_t>(v & 0xFF);
            v >>= 8;
        }
    } else {
        for (int i = 0; i < 8; ++i) {
            q[i] = static_cast<uint8_t>(v & 0xFF);
            v >>= 8;
        }
    }
}

typedef struct sha3_context_ {
    uint64_t s[SHA3_KECCAK_SPONGE_WORDS];
    uint64_t saved;
    uint32_t byteIndex;
    uint32_t wordIndex;
    uint32_t capacityWords;
} sha3_context;

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

static void keccakf(uint64_t s[25]) {
    int i, j, round;
    uint64_t t, bc[5];

    for (round = 0; round < SHA3_KECCAK_ROUNDS; round++) {
        // Theta
        for (i = 0; i < 5; i++) {
            bc[i] = s[i] ^ s[i + 5] ^ s[i + 10] ^ s[i + 15] ^ s[i + 20];
        }
        for (i = 0; i < 5; i++) {
            t = bc[(i + 4) % 5] ^ ROTL64(bc[(i + 1) % 5], 1);
            for (j = 0; j < 25; j += 5) {
                s[j + i] ^= t;
            }
        }

        // Rho and Pi
        t = s[1];
        for (i = 0; i < 24; i++) {
            j = keccakf_piln[i];
            bc[0] = s[j];
            s[j] = ROTL64(t, keccakf_rotc[i]);
            t = bc[0];
        }

        // Chi
        for (j = 0; j < 25; j += 5) {
            for (i = 0; i < 5; i++) {
                bc[i] = s[j + i];
            }
            for (i = 0; i < 5; i++) {
                s[j + i] ^= (~bc[(i + 1) % 5]) & bc[(i + 2) % 5];
            }
        }

        // Iota
        s[0] ^= keccakf_rndc[round];
    }
}

static void sha3_Init(sha3_context* ctx, unsigned bitSize) {
    assert(bitSize == 256 || bitSize == 384 || bitSize == 512);
    std::memset(ctx, 0, sizeof(*ctx));
    ctx->capacityWords = 2 * bitSize / (8 * sizeof(uint64_t));
}

static void sha3_Seed(sha3_context* ctx, uint64_t seed) {
    if (ctx->capacityWords >= 2) {
        ctx->s[SHA3_KECCAK_SPONGE_WORDS - 2] ^= seed;
        ctx->s[SHA3_KECCAK_SPONGE_WORDS - 1] ^= seed * UINT64_C(0x9E3779B97F4A7C15);
    } else {
        ctx->s[SHA3_KECCAK_SPONGE_WORDS - 1] ^= seed;
    }
}

template <bool bswap>
static void sha3_Process(sha3_context* ctx, const uint8_t* in, size_t inlen) {
    uint32_t old_tail = (8 - ctx->byteIndex) & 7;
    uint32_t tail;
    size_t words, i;

    if (inlen == 0) return;

    if (inlen < old_tail) {
        while (inlen--) {
            ctx->saved |= static_cast<uint64_t>(*(in++)) << ((ctx->byteIndex++) * 8);
        }
        return;
    }

    if (old_tail) {
        inlen -= old_tail;
        while (old_tail--) {
            ctx->saved |= static_cast<uint64_t>(*(in++)) << ((ctx->byteIndex++) * 8);
        }
        ctx->s[ctx->wordIndex] ^= ctx->saved;
        ctx->byteIndex = 0;
        ctx->saved = 0;

        if (++ctx->wordIndex == (SHA3_KECCAK_SPONGE_WORDS - ctx->capacityWords)) {
            keccakf(ctx->s);
            ctx->wordIndex = 0;
        }
    }

    words = inlen / sizeof(uint64_t);
    tail = static_cast<uint32_t>(inlen - words * sizeof(uint64_t));

    for (i = 0; i < words; i++, in += sizeof(uint64_t)) {
        uint64_t t = GET_U64<bswap>(in, 0);
        ctx->s[ctx->wordIndex] ^= t;

        if (++ctx->wordIndex == (SHA3_KECCAK_SPONGE_WORDS - ctx->capacityWords)) {
            keccakf(ctx->s);
            ctx->wordIndex = 0;
        }
    }

    while (tail--) {
        ctx->saved |= static_cast<uint64_t>(*(in++)) << ((ctx->byteIndex++) * 8);
    }
}

template <bool bswap>
static void sha3_Finalize(sha3_context* ctx, size_t digest_words, uint8_t* digest) {
    // SHA-3 domain suffix 01 plus padding, encoded as 0x06.
    uint64_t t = static_cast<uint64_t>(
        (static_cast<uint64_t>(0x02 | (1 << 2))) << ((ctx->byteIndex) * 8)
    );

    ctx->s[ctx->wordIndex] ^= ctx->saved ^ t;
    ctx->s[SHA3_KECCAK_SPONGE_WORDS - ctx->capacityWords - 1] ^= UINT64_C(0x8000000000000000);

    keccakf(ctx->s);

    uint32_t maxdigest_words = ctx->capacityWords / 2;
    if (digest_words > maxdigest_words) {
        digest_words = maxdigest_words;
    }

    for (size_t i = 0; i < digest_words; i++) {
        PUT_U64<bswap>(ctx->s[i], digest, 8 * i);
    }
}

template <uint32_t hashbits, bool bswap>
static void SHA3_256_core(const void* in, const size_t len, const seed_t seed, void* out) {
    sha3_context context;
    sha3_Init(&context, 256);
    sha3_Seed(&context, static_cast<uint64_t>(seed));
    sha3_Process<bswap>(&context, static_cast<const uint8_t*>(in), len);
    sha3_Finalize<bswap>(&context, (hashbits + 63) / 64, static_cast<uint8_t*>(out));
}

static void hash_sha3_256_bytes(const uint8_t* message, size_t len, uint8_t digest[32]) {
    // seed = 0 preserves unseeded SHA-3 behavior, matching the prior experiment setup.
    SHA3_256_core<256, false>(message, len, 0, digest);
}

// ============================================================
// Generator utilities
// ============================================================
static uint64_t parse_u64(const std::string& s) {
    std::size_t pos = 0;
    uint64_t v = std::stoull(s, &pos, 0);
    if (pos != s.size()) throw std::runtime_error("invalid integer: " + s);
    return v;
}

static void put_le64(uint64_t x, uint8_t* out) {
    for (int i = 0; i < 8; ++i) {
        out[i] = static_cast<uint8_t>((x >> (8 * i)) & 0xFFU);
    }
}

static uint64_t splitmix64_permute(uint64_t x) {
    x += UINT64_C(0x9E3779B97F4A7C15);
    x = (x ^ (x >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
    x = (x ^ (x >> 27)) * UINT64_C(0x94D049BB133111EB);
    x = x ^ (x >> 31);
    return x;
}

enum class OutputFormat {
    Binary,
    Ascii
};

struct Options {
    uint64_t sequences = 1000;
    uint64_t bits_per_seq = 1000000;
    OutputFormat format = OutputFormat::Binary;
    std::string output = "sha3_256_nist_1000x1M.bin";
    uint64_t domain_salt = UINT64_C(0x534841335f4e4953); // ASCII-ish: SHA3_NIS
    uint64_t progress_every = 10;
};

static void print_usage(const char* prog) {
    std::cerr
        << "Usage: " << prog << " [options]\n\n"
        << "Options:\n"
        << "  --sequences N          Number of NIST STS bitstreams, default 1000\n"
        << "  --bits-per-seq N       Bits per bitstream, default 1000000\n"
        << "  --format binary|ascii  Output format, default binary\n"
        << "  --output FILE          Output file path\n"
        << "  --domain-salt N        64-bit deterministic domain salt\n"
        << "  --progress-every N     Print progress every N sequences, default 10\n"
        << "  --help                 Show this help\n\n"
        << "Default full dataset: 1000 x 1,000,000 bits = 1,000,000,000 bits.\n";
}

static Options parse_args(int argc, char** argv) {
    Options opt;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto require_value = [&](const std::string& flag) -> std::string {
            if (i + 1 >= argc) throw std::runtime_error("missing value for " + flag);
            return argv[++i];
        };

        if (a == "--sequences") {
            opt.sequences = parse_u64(require_value(a));
        } else if (a == "--bits-per-seq") {
            opt.bits_per_seq = parse_u64(require_value(a));
        } else if (a == "--format") {
            std::string v = require_value(a);
            if (v == "binary" || v == "bin") opt.format = OutputFormat::Binary;
            else if (v == "ascii" || v == "txt") opt.format = OutputFormat::Ascii;
            else throw std::runtime_error("invalid --format: " + v + " (use binary or ascii)");
        } else if (a == "--output") {
            opt.output = require_value(a);
        } else if (a == "--domain-salt") {
            opt.domain_salt = parse_u64(require_value(a));
        } else if (a == "--progress-every") {
            opt.progress_every = parse_u64(require_value(a));
        } else if (a == "--help" || a == "-h") {
            print_usage(argv[0]);
            std::exit(0);
        } else {
            throw std::runtime_error("unknown option: " + a);
        }
    }

    if (opt.sequences == 0) throw std::runtime_error("--sequences must be > 0");
    if (opt.bits_per_seq == 0) throw std::runtime_error("--bits-per-seq must be > 0");
    return opt;
}

class BitWriter {
public:
    BitWriter(std::ofstream& out, OutputFormat format)
        : out_(out), format_(format) {}

    void write_bit(int bit) {
        bit &= 1;
        if (format_ == OutputFormat::Ascii) {
            char c = bit ? '1' : '0';
            out_.write(&c, 1);
            ++bits_written_;
            return;
        }

        current_byte_ = static_cast<uint8_t>((current_byte_ << 1) | static_cast<uint8_t>(bit));
        ++bits_in_current_;
        ++bits_written_;
        if (bits_in_current_ == 8) {
            out_.put(static_cast<char>(current_byte_));
            current_byte_ = 0;
            bits_in_current_ = 0;
        }
    }

    void write_digest_bits(const uint8_t digest[32], uint64_t max_bits) {
        uint64_t written = 0;
        for (int i = 0; i < 32 && written < max_bits; ++i) {
            for (int b = 7; b >= 0 && written < max_bits; --b) {
                write_bit((digest[i] >> b) & 1U);
                ++written;
            }
        }
    }

    void flush() {
        if (format_ == OutputFormat::Binary && bits_in_current_ != 0) {
            current_byte_ <<= (8 - bits_in_current_);
            out_.put(static_cast<char>(current_byte_));
            current_byte_ = 0;
            bits_in_current_ = 0;
        }
    }

    uint64_t bits_written() const { return bits_written_; }

private:
    std::ofstream& out_;
    OutputFormat format_;
    uint8_t current_byte_ = 0;
    int bits_in_current_ = 0;
    uint64_t bits_written_ = 0;
};

static void build_message(uint64_t domain_salt, uint64_t sequence_index, uint64_t block_counter, uint8_t msg[24]) {
    const uint64_t seq_salt = splitmix64_permute(domain_salt ^ (UINT64_C(0x9E3779B97F4A7C15) * (sequence_index + 1)));
    put_le64(seq_salt, msg + 0);
    put_le64(sequence_index, msg + 8);
    put_le64(block_counter, msg + 16);
}

int main(int argc, char** argv) {
    try {
        Options opt = parse_args(argc, argv);

        std::ofstream out(opt.output, std::ios::binary);
        if (!out) throw std::runtime_error("failed to open output file: " + opt.output);

        BitWriter writer(out, opt.format);
        const auto start = Clock::now();

        std::cout << "============================================================\n";
        std::cout << "  SHA-3-256 NIST STS BITSTREAM GENERATOR\n";
        std::cout << "============================================================\n";
        std::cout << "sequences        : " << opt.sequences << "\n";
        std::cout << "bits per sequence: " << opt.bits_per_seq << "\n";
        std::cout << "total bits       : " << (opt.sequences * opt.bits_per_seq) << "\n";
        std::cout << "format           : " << (opt.format == OutputFormat::Binary ? "binary" : "ascii") << "\n";
        std::cout << "output           : " << opt.output << "\n";
        std::cout << "domain salt      : 0x" << std::hex << opt.domain_salt << std::dec << "\n\n";

        uint8_t msg[24];
        uint8_t digest[32];

        for (uint64_t seq = 0; seq < opt.sequences; ++seq) {
            uint64_t produced = 0;
            uint64_t block = 0;

            while (produced < opt.bits_per_seq) {
                build_message(opt.domain_salt, seq, block, msg);
                hash_sha3_256_bytes(msg, sizeof(msg), digest);

                const uint64_t remaining = opt.bits_per_seq - produced;
                const uint64_t take = remaining < 256 ? remaining : 256;
                writer.write_digest_bits(digest, take);
                produced += take;
                ++block;
            }

            if (opt.progress_every > 0 && ((seq + 1) % opt.progress_every == 0 || seq + 1 == opt.sequences)) {
                const auto now = Clock::now();
                const double seconds = std::chrono::duration<double>(now - start).count();
                std::cout << "[progress] sequence " << (seq + 1) << "/" << opt.sequences
                          << " | bits=" << writer.bits_written()
                          << " | seconds=" << std::fixed << std::setprecision(2) << seconds << "\n";
            }
        }

        writer.flush();
        out.close();

        const auto end = Clock::now();
        const double seconds = std::chrono::duration<double>(end - start).count();

        std::cout << "\nDone.\n";
        std::cout << "written bits : " << writer.bits_written() << "\n";
        std::cout << "output file  : " << opt.output << "\n";
        std::cout << "time seconds : " << std::fixed << std::setprecision(3) << seconds << "\n";

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
