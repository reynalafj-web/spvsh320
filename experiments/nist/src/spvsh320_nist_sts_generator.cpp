#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// ============================================================================
//  NIST STS INPUT GENERATOR: SpVSH-320-FR2
// ============================================================================
//
//  Purpose:
//    Generate binary or ASCII bitstreams for NIST SP 800-22/STSuite testing.
//
//  Hash construction:
//    SpVSH-320-FR2, using the finalized SpVSH-320 algorithm with:
//      - 320-bit internal state: five 64-bit lanes
//      - rate lane S0
//      - original PermuteSpVSH transformation
//      - 0x80 padding block
//      - two output-finalization rounds by default before first output word
//      - 256-bit digest generated from four S0 squeeze words
//
//  Message schedule for bitstream generation:
//    M_{seq,ctr} = LE64(domain_salt) || LE64(seq_index) || LE64(block_counter)
//
//  The generator concatenates digest bytes until each NIST bitstream contains
//  exactly bits_per_sequence bits. For binary output, bits_per_sequence must be
//  a multiple of 8.
// ============================================================================

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

static inline uint64_t MulModP1(uint64_t a, uint64_t b) { return MulModP(a, b, P1, C1); }
static inline uint64_t MulModP2(uint64_t a, uint64_t b) { return MulModP(a, b, P2, C2); }

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

    state.s[0] = t0 ^ (t1 >> 17) ^ (t1 << 47);
    state.s[1] = t1 ^ (t2 >> 11) ^ (t2 << 53);
    state.s[2] = t2 ^ (t3 >> 19) ^ (t3 << 45);
    state.s[3] = t3 ^ (t4 >> 13) ^ (t4 << 51);
    state.s[4] = t4 ^ (t0 >> 23) ^ (t0 << 41);
}

static inline uint64_t GET_U64_LE(const uint8_t* ptr) {
    uint64_t val;
    std::memcpy(&val, ptr, 8);
    return val;
}

static inline void PUT_U64_LE(uint64_t v, uint8_t* out) {
    for (int i = 0; i < 8; ++i) {
        out[i] = static_cast<uint8_t>((v >> (8 * i)) & 0xFFU);
    }
}

static inline void write_le64(uint64_t v, uint8_t* out) {
    PUT_U64_LE(v, out);
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

static inline void ApplyFinalizationRounds(PentaState& state, unsigned rounds) {
    for (unsigned r = 0; r < rounds; ++r) {
        PermuteSpVSH(state);
    }
}

static void SpVSH320_FR2_256(const uint8_t* in, size_t len, uint8_t* out, unsigned finalization_rounds) {
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
        std::memcpy(buffer, in, length);
    }
    buffer[length] = 0x80;
    state.s[0] ^= GET_U64_LE(buffer);
    PermuteSpVSH(state);

    ApplyFinalizationRounds(state, finalization_rounds);

    uint64_t h[4];
    for (int i = 0; i < 4; ++i) {
        h[i] = state.s[0];
        if (i < 3) {
            PermuteSpVSH(state);
        }
    }

    // Canonical little-endian 64-bit word serialization, matching the original
    // implementation's memcpy behavior on x86-64 while remaining explicit.
    for (int i = 0; i < 4; ++i) {
        PUT_U64_LE(h[i], out + 8 * i);
    }
}

struct Options {
    uint64_t sequences = 1000;
    uint64_t bits_per_sequence = 1000000;
    uint64_t start_sequence = 0;
    uint64_t domain_salt = UINT64_C(0x5350565348333230); // ASCII-like SPVSH320
    unsigned finalization_rounds = 2;
    uint64_t progress_every = 10;
    std::string output_path = "spvsh320_fr2_nist_1000x1M.bin";
    std::string format = "binary"; // binary or ascii
};

static uint64_t parse_u64(const std::string& s) {
    std::size_t pos = 0;
    uint64_t v = std::stoull(s, &pos, 0);
    if (pos != s.size()) throw std::runtime_error("invalid integer: " + s);
    return v;
}

static void print_usage(const char* prog) {
    std::cerr
        << "Usage: " << prog << " [options]\n\n"
        << "Options:\n"
        << "  --sequences N              Number of NIST bitstreams, default 1000\n"
        << "  --bits-per-seq N           Bits per sequence, default 1000000\n"
        << "  --start-sequence N         First sequence index, default 0\n"
        << "  --domain-salt N            64-bit domain salt, default 0x5350565348333230\n"
        << "  --finalization-rounds N    Extra finalization rounds, default 2\n"
        << "  --format binary|ascii      Output format for STS, default binary\n"
        << "  --output PATH              Output bitstream path\n"
        << "  --progress-every N         Progress interval in sequences, default 10\n"
        << "  --help                     Show this help\n\n"
        << "Example:\n"
        << "  " << prog << " --sequences 1000 --bits-per-seq 1000000 --format binary --output data/spvsh320_fr2_1000x1M.bin\n";
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
        } else if (a == "--bits-per-seq" || a == "--bits-per-sequence") {
            opt.bits_per_sequence = parse_u64(require_value(a));
        } else if (a == "--start-sequence") {
            opt.start_sequence = parse_u64(require_value(a));
        } else if (a == "--domain-salt") {
            opt.domain_salt = parse_u64(require_value(a));
        } else if (a == "--finalization-rounds" || a == "--final-rounds") {
            opt.finalization_rounds = static_cast<unsigned>(std::stoul(require_value(a)));
        } else if (a == "--format") {
            opt.format = require_value(a);
            if (opt.format != "binary" && opt.format != "ascii") {
                throw std::runtime_error("--format must be binary or ascii");
            }
        } else if (a == "--output") {
            opt.output_path = require_value(a);
        } else if (a == "--progress-every") {
            opt.progress_every = parse_u64(require_value(a));
        } else if (a == "--help" || a == "-h") {
            print_usage(argv[0]);
            std::exit(0);
        } else {
            throw std::runtime_error("unknown option: " + a);
        }
    }

    if (opt.sequences == 0) throw std::runtime_error("--sequences must be positive");
    if (opt.bits_per_sequence == 0) throw std::runtime_error("--bits-per-seq must be positive");
    if (opt.format == "binary" && (opt.bits_per_sequence % 8 != 0)) {
        throw std::runtime_error("binary output requires bits_per_sequence to be a multiple of 8");
    }
    return opt;
}

static void write_ascii_bits(std::ofstream& out, const uint8_t* bytes, size_t byte_count, uint64_t& bits_remaining) {
    for (size_t i = 0; i < byte_count && bits_remaining > 0; ++i) {
        for (int bit = 7; bit >= 0 && bits_remaining > 0; --bit) {
            out.put(((bytes[i] >> bit) & 1U) ? '1' : '0');
            --bits_remaining;
        }
    }
}

static void write_metadata(const Options& opt, uint64_t total_bits, uint64_t total_hashes, double seconds) {
    std::ofstream meta(opt.output_path + ".meta.txt");
    if (!meta) return;
    meta << "NIST STS Input Generator: SpVSH-320-FR2\n";
    meta << "output_path=" << opt.output_path << "\n";
    meta << "format=" << opt.format << "\n";
    meta << "sequences=" << opt.sequences << "\n";
    meta << "bits_per_sequence=" << opt.bits_per_sequence << "\n";
    meta << "total_bits=" << total_bits << "\n";
    meta << "total_bytes=" << ((total_bits + 7) / 8) << "\n";
    meta << "total_hashes=" << total_hashes << "\n";
    meta << "message_format=LE64(domain_salt)||LE64(sequence_index)||LE64(block_counter)\n";
    meta << "domain_salt=0x" << std::hex << opt.domain_salt << std::dec << "\n";
    meta << "start_sequence=" << opt.start_sequence << "\n";
    meta << "digest_bits=256\n";
    meta << "finalization_rounds=" << opt.finalization_rounds << "\n";
    meta << "seconds=" << std::fixed << std::setprecision(6) << seconds << "\n";
}

int main(int argc, char** argv) {
    try {
        Options opt = parse_args(argc, argv);
        std::ofstream out(opt.output_path, std::ios::binary);
        if (!out) throw std::runtime_error("failed to open output file: " + opt.output_path);

        const auto start = std::chrono::high_resolution_clock::now();
        const uint64_t total_bits = opt.sequences * opt.bits_per_sequence;
        uint64_t total_hashes = 0;

        std::cout << "============================================================\n";
        std::cout << "  NIST STS INPUT GENERATOR: SpVSH-320-FR2\n";
        std::cout << "============================================================\n";
        std::cout << "Output              : " << opt.output_path << "\n";
        std::cout << "Format              : " << opt.format << "\n";
        std::cout << "Sequences           : " << opt.sequences << "\n";
        std::cout << "Bits per sequence   : " << opt.bits_per_sequence << "\n";
        std::cout << "Total bits          : " << total_bits << "\n";
        std::cout << "Finalization rounds : " << opt.finalization_rounds << "\n";
        std::cout << "Domain salt         : 0x" << std::hex << opt.domain_salt << std::dec << "\n\n";

        uint8_t message[24];
        uint8_t digest[32];

        for (uint64_t seq = 0; seq < opt.sequences; ++seq) {
            uint64_t bits_remaining = opt.bits_per_sequence;
            uint64_t bytes_remaining = opt.format == "binary" ? (opt.bits_per_sequence / 8) : 0;
            uint64_t block_counter = 0;
            const uint64_t seq_index = opt.start_sequence + seq;

            while (bits_remaining > 0) {
                write_le64(opt.domain_salt, message + 0);
                write_le64(seq_index, message + 8);
                write_le64(block_counter, message + 16);

                SpVSH320_FR2_256(message, sizeof(message), digest, opt.finalization_rounds);
                ++total_hashes;
                ++block_counter;

                if (opt.format == "binary") {
                    const size_t take = static_cast<size_t>(std::min<uint64_t>(32, bytes_remaining));
                    out.write(reinterpret_cast<const char*>(digest), static_cast<std::streamsize>(take));
                    bytes_remaining -= take;
                    bits_remaining -= static_cast<uint64_t>(take) * 8ULL;
                } else {
                    write_ascii_bits(out, digest, sizeof(digest), bits_remaining);
                }
            }

            if (opt.progress_every != 0 && ((seq + 1) % opt.progress_every == 0 || seq + 1 == opt.sequences)) {
                std::cout << "Generated sequence " << (seq + 1) << "/" << opt.sequences << "\n";
            }
        }

        out.close();
        const auto end = std::chrono::high_resolution_clock::now();
        const double seconds = std::chrono::duration<double>(end - start).count();
        write_metadata(opt, total_bits, total_hashes, seconds);

        std::cout << "\nDone.\n";
        std::cout << "Total hashes         : " << total_hashes << "\n";
        std::cout << "Elapsed seconds      : " << std::fixed << std::setprecision(4) << seconds << "\n";
        std::cout << "Metadata             : " << opt.output_path << ".meta.txt\n";
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
