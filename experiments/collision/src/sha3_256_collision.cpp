#include <algorithm>
#include <atomic>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

using Clock = std::chrono::high_resolution_clock;
using seed_t = uint64_t;

// ============================================================================
//  MULTICORE TRUNCATED COLLISION EXPERIMENT: SHA-3-256
// ============================================================================
//
//  Scenario aligned with the revised VSH-256 and SpVSH-320 collision experiment:
//
//    F_t(x) = rightmost_t(SHA3-256(LE64(salt_{t,r}) || LE64(x)))
//
//  Collision search method:
//    Pollard's rho with Floyd's tortoise-and-hare cycle finding.
//
//  Default experiment:
//    truncation bits : 16,24,32,40,48,56
//    repetitions     : 50 for 16/24/32-bit, 30 for 40-bit,
//                      10 for 48-bit, 5 for 56-bit
//    phase-1 limit   : 2^32 by default
//
//  The SHA-3-256 core below is retained from the preimage experiment. Only the
//  experimental framework is changed from truncated preimage to truncated
//  collision. The hash algorithm itself is not modified.
// ============================================================================

// ============================================================
// SHA3-256 core implementation
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
        for (int i = 0; i < 8; ++i) v = (v << 8) | q[i];
    } else {
        for (int i = 7; i >= 0; --i) v = (v << 8) | q[i];
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
        for (i = 0; i < 5; i++) {
            bc[i] = s[i] ^ s[i + 5] ^ s[i + 10] ^ s[i + 15] ^ s[i + 20];
        }

        for (i = 0; i < 5; i++) {
            t = bc[(i + 4) % 5] ^ ROTL64(bc[(i + 1) % 5], 1);
            for (j = 0; j < 25; j += 5) s[j + i] ^= t;
        }

        t = s[1];
        for (i = 0; i < 24; i++) {
            j = keccakf_piln[i];
            bc[0] = s[j];
            s[j] = ROTL64(t, keccakf_rotc[i]);
            t = bc[0];
        }

        for (j = 0; j < 25; j += 5) {
            for (i = 0; i < 5; i++) bc[i] = s[j + i];
            for (i = 0; i < 5; i++) s[j + i] ^= (~bc[(i + 1) % 5]) & bc[(i + 2) % 5];
        }

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
        while (inlen--) ctx->saved |= static_cast<uint64_t>(*(in++)) << ((ctx->byteIndex++) * 8);
        return;
    }

    if (old_tail) {
        inlen -= old_tail;
        while (old_tail--) ctx->saved |= static_cast<uint64_t>(*(in++)) << ((ctx->byteIndex++) * 8);
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
        uint64_t tt = GET_U64<bswap>(in, 0);
        ctx->s[ctx->wordIndex] ^= tt;
        if (++ctx->wordIndex == (SHA3_KECCAK_SPONGE_WORDS - ctx->capacityWords)) {
            keccakf(ctx->s);
            ctx->wordIndex = 0;
        }
    }

    while (tail--) ctx->saved |= static_cast<uint64_t>(*(in++)) << ((ctx->byteIndex++) * 8);
}

template <bool bswap>
static void sha3_Finalize(sha3_context* ctx, size_t digest_words, uint8_t* digest) {
    uint64_t t = static_cast<uint64_t>((static_cast<uint64_t>(0x02 | (1 << 2))) << ((ctx->byteIndex) * 8));
    ctx->s[ctx->wordIndex] ^= ctx->saved ^ t;
    ctx->s[SHA3_KECCAK_SPONGE_WORDS - ctx->capacityWords - 1] ^= UINT64_C(0x8000000000000000);

    keccakf(ctx->s);

    uint32_t maxdigest_words = ctx->capacityWords / 2;
    if (digest_words > maxdigest_words) digest_words = maxdigest_words;

    for (size_t i = 0; i < digest_words; i++) PUT_U64<bswap>(ctx->s[i], digest, 8 * i);
}

template <uint32_t hashbits, bool bswap>
static void SHA3_256_core(const void* in, const size_t len, const seed_t seed, void* out) {
    sha3_context context;
    sha3_Init(&context, 256);
    sha3_Seed(&context, static_cast<uint64_t>(seed));
    sha3_Process<bswap>(&context, static_cast<const uint8_t*>(in), len);
    sha3_Finalize<bswap>(&context, (hashbits + 63) / 64, static_cast<uint8_t*>(out));
}

// ============================================================================
// Experiment utilities
// ============================================================================
static inline uint64_t mask_bits(unsigned t) {
    if (t == 64) return UINT64_MAX;
    return (UINT64_C(1) << t) - 1;
}

static inline uint64_t rightmost_bits_from_digest(const uint8_t digest[32], unsigned t) {
    if (t == 0 || t > 64) throw std::runtime_error("truncation bits must be in 1..64");
    uint64_t last64 = 0;
    for (int i = 24; i < 32; ++i) last64 = (last64 << 8) | digest[i];
    return t == 64 ? last64 : (last64 & mask_bits(t));
}

static inline void put_u64_le(uint8_t* out, uint64_t x) {
    for (int i = 0; i < 8; ++i) out[i] = static_cast<uint8_t>((x >> (8 * i)) & 0xFFU);
}

static std::vector<uint8_t> collision_message_bytes(uint64_t salt, uint64_t state) {
    std::vector<uint8_t> out(16, 0);
    put_u64_le(out.data(), salt);
    put_u64_le(out.data() + 8, state);
    return out;
}

static std::string bytes_to_hex(const std::vector<uint8_t>& bytes) {
    std::ostringstream oss;
    oss << std::hex << std::setfill('0');
    for (uint8_t b : bytes) oss << std::setw(2) << static_cast<unsigned>(b);
    return oss.str();
}

static std::string truncated_hex(uint64_t value, unsigned t) {
    const unsigned width = (t + 3) / 4;
    std::ostringstream oss;
    oss << "0x" << std::hex << std::setfill('0') << std::setw(width) << value;
    return oss.str();
}

static std::string u64_hex(uint64_t value, unsigned width = 16) {
    std::ostringstream oss;
    oss << "0x" << std::hex << std::setfill('0') << std::setw(width) << value;
    return oss.str();
}

static std::vector<unsigned> parse_bits_list(const std::string& s) {
    std::vector<unsigned> out;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, ',')) {
        if (item.empty()) continue;
        unsigned v = static_cast<unsigned>(std::stoul(item));
        if (v == 0 || v > 64) throw std::runtime_error("truncation bits must be in the range 1..64");
        out.push_back(v);
    }
    if (out.empty()) throw std::runtime_error("empty --bits list");
    return out;
}

static uint64_t parse_u64(const std::string& s) {
    std::size_t pos = 0;
    uint64_t v = std::stoull(s, &pos, 0);
    if (pos != s.size()) throw std::runtime_error("invalid integer: " + s);
    return v;
}

static inline uint64_t splitmix64(uint64_t x) {
    x += UINT64_C(0x9E3779B97F4A7C15);
    x = (x ^ (x >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
    x = (x ^ (x >> 27)) * UINT64_C(0x94D049BB133111EB);
    return x ^ (x >> 31);
}

static uint64_t salt_for(unsigned t, int repetition) {
    uint64_t x = UINT64_C(0xC0111510D15EA5E5);
    x ^= (static_cast<uint64_t>(t) << 48);
    x ^= (static_cast<uint64_t>(repetition) * UINT64_C(0x9E3779B97F4A7C15));
    return splitmix64(x);
}

static uint64_t seed_state_for(unsigned t, int repetition, int seed_attempt) {
    uint64_t x = UINT64_C(0xA0761D6478BD642F);
    x ^= (static_cast<uint64_t>(t) << 40);
    x ^= (static_cast<uint64_t>(repetition) << 16);
    x ^= static_cast<uint64_t>(seed_attempt);
    return splitmix64(x) & mask_bits(t);
}

static uint64_t transition_sha3_rightproj(uint64_t state, uint64_t salt, unsigned t, uint64_t& hash_evaluations) {
    uint8_t msg[16] = {0};
    put_u64_le(msg, salt);
    put_u64_le(msg + 8, state);

    uint8_t digest[32] = {0};
    SHA3_256_core<256, false>(msg, sizeof(msg), 0, digest);
    ++hash_evaluations;
    return rightmost_bits_from_digest(digest, t);
}

static long double expected_collision_steps(unsigned t) {
    static const long double sqrt_pi_over_2 = 1.2533141373155002512L;
    return sqrt_pi_over_2 * std::pow(2.0L, static_cast<long double>(t) / 2.0L);
}

static int repetitions_for_truncation(unsigned t) {
    if (t <= 32) return 50;
    if (t == 40) return 30;
    if (t == 48) return 10;
    if (t == 56) return 5;
    return 3;
}

struct TrialResult {
    unsigned trunc_bits = 0;
    int repetition = 0;
    int seed_attempt = 0;
    std::string seed_state_hex;
    bool found = false;
    uint64_t rho_meet_steps = 0;
    uint64_t mu = 0;
    uint64_t lambda = 0;
    uint64_t collision_index = 0;
    uint64_t hash_evaluations = 0;
    long double expected_steps = 0.0L;
    double log2_collision_index = 0.0;
    long double ratio_to_expected = 0.0L;
    double seconds = 0.0;
    std::string input_a_hex;
    std::string input_b_hex;
    std::string collision_trunc_hex;
};

struct SummaryResult {
    unsigned trunc_bits = 0;
    int repetitions = 0;
    int success_count = 0;
    long double expected_steps = 0.0L;
    long double mean_collision_index = 0.0L;
    long double median_collision_index = 0.0L;
    uint64_t min_collision_index = 0;
    uint64_t max_collision_index = 0;
    long double stddev_collision_index = 0.0L;
    long double mean_hash_evaluations = 0.0L;
    long double mean_log2_collision_index = 0.0L;
    long double mean_ratio_to_expected = 0.0L;
    long double mean_seconds = 0.0L;
};

struct Options {
    std::vector<unsigned> bits = {16, 24, 32, 40, 48, 56};
    unsigned workers = 0;
    uint64_t max_phase_steps = (uint64_t(1) << 32);
    int max_seed_attempts = 16;
    std::string out_prefix = "sha3_256_collision_rightproj";
};

static void print_usage(const char* prog) {
    std::cerr
        << "Usage: " << prog << " [options]\n\n"
        << "Options:\n"
        << "  --bits 16,24,32,40,48,56  Comma-separated truncation levels\n"
        << "  --workers N                 Parallel repetition workers; 0 = auto\n"
        << "  --max-phase N               Maximum Floyd phase-1 meeting iterations\n"
        << "  --max-phase-power P         Set max phase iterations to 2^P\n"
        << "  --max-seed-attempts N       Deterministic seed retries when mu=0 or verification fails\n"
        << "  --out-prefix PREFIX         Output filename prefix\n"
        << "  --help                      Show this help\n\n"
        << "Example:\n"
        << "  " << prog << " --workers 50 --bits 16,24,32,40,48,56 --out-prefix sha3_collision_run01\n";
}

static Options parse_args(int argc, char** argv) {
    Options opt;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto require_value = [&](const std::string& flag) -> std::string {
            if (i + 1 >= argc) throw std::runtime_error("missing value for " + flag);
            return argv[++i];
        };
        if (a == "--bits") {
            opt.bits = parse_bits_list(require_value(a));
        } else if (a == "--workers") {
            opt.workers = static_cast<unsigned>(std::stoul(require_value(a)));
        } else if (a == "--max-phase") {
            opt.max_phase_steps = parse_u64(require_value(a));
        } else if (a == "--max-phase-power") {
            unsigned p = static_cast<unsigned>(std::stoul(require_value(a)));
            if (p >= 63) throw std::runtime_error("--max-phase-power must be < 63");
            opt.max_phase_steps = (uint64_t(1) << p);
        } else if (a == "--max-seed-attempts") {
            opt.max_seed_attempts = std::stoi(require_value(a));
            if (opt.max_seed_attempts <= 0) throw std::runtime_error("--max-seed-attempts must be positive");
        } else if (a == "--out-prefix") {
            opt.out_prefix = require_value(a);
        } else if (a == "--help" || a == "-h") {
            print_usage(argv[0]);
            std::exit(0);
        } else {
            throw std::runtime_error("unknown option: " + a);
        }
    }
    if (opt.workers == 0) {
        opt.workers = std::thread::hardware_concurrency();
        if (opt.workers == 0) opt.workers = 1;
    }
    return opt;
}

static TrialResult run_one_trial(unsigned t, int repetition, const Options& opt) {
    const uint64_t salt = salt_for(t, repetition);
    const long double expected = expected_collision_steps(t);

    for (int seed_attempt = 0; seed_attempt < opt.max_seed_attempts; ++seed_attempt) {
        uint64_t evals = 0;
        const uint64_t x0 = seed_state_for(t, repetition, seed_attempt);
        auto start = Clock::now();

        auto F = [&](uint64_t x) -> uint64_t {
            return transition_sha3_rightproj(x, salt, t, evals);
        };

        uint64_t tortoise = F(x0);
        uint64_t hare = F(F(x0));
        uint64_t meet_steps = 1;

        while (tortoise != hare && meet_steps < opt.max_phase_steps) {
            tortoise = F(tortoise);
            hare = F(F(hare));
            ++meet_steps;
        }

        if (tortoise != hare) {
            auto end = Clock::now();
            TrialResult r;
            r.trunc_bits = t;
            r.repetition = repetition;
            r.seed_attempt = seed_attempt;
            r.seed_state_hex = u64_hex(x0, (t + 3) / 4);
            r.found = false;
            r.rho_meet_steps = meet_steps;
            r.hash_evaluations = evals;
            r.expected_steps = expected;
            r.seconds = std::chrono::duration<double>(end - start).count();
            return r;
        }

        uint64_t mu = 0;
        tortoise = x0;
        while (tortoise != hare) {
            tortoise = F(tortoise);
            hare = F(hare);
            ++mu;
        }

        uint64_t lambda = 1;
        hare = F(tortoise);
        while (tortoise != hare) {
            hare = F(hare);
            ++lambda;
        }

        if (mu == 0) {
            continue;
        }

        auto state_at = [&](uint64_t idx) -> uint64_t {
            uint64_t x = x0;
            for (uint64_t i = 0; i < idx; ++i) x = F(x);
            return x;
        };

        const uint64_t xa = state_at(mu - 1);
        const uint64_t xb = state_at(mu + lambda - 1);
        const uint64_t ya = F(xa);
        const uint64_t yb = F(xb);

        if (xa == xb || ya != yb) {
            continue;
        }

        auto end = Clock::now();
        TrialResult r;
        r.trunc_bits = t;
        r.repetition = repetition;
        r.seed_attempt = seed_attempt;
        r.seed_state_hex = u64_hex(x0, (t + 3) / 4);
        r.found = true;
        r.rho_meet_steps = meet_steps;
        r.mu = mu;
        r.lambda = lambda;
        r.collision_index = mu + lambda;
        r.hash_evaluations = evals;
        r.expected_steps = expected;
        r.log2_collision_index = std::log2(static_cast<double>(r.collision_index));
        r.ratio_to_expected = static_cast<long double>(r.collision_index) / expected;
        r.seconds = std::chrono::duration<double>(end - start).count();
        r.input_a_hex = bytes_to_hex(collision_message_bytes(salt, xa));
        r.input_b_hex = bytes_to_hex(collision_message_bytes(salt, xb));
        r.collision_trunc_hex = truncated_hex(ya, t);
        return r;
    }

    TrialResult r;
    r.trunc_bits = t;
    r.repetition = repetition;
    r.seed_attempt = opt.max_seed_attempts;
    r.found = false;
    r.expected_steps = expected;
    return r;
}

static SummaryResult summarize(unsigned t, const std::vector<TrialResult>& all_results) {
    SummaryResult s;
    s.trunc_bits = t;
    s.expected_steps = expected_collision_steps(t);

    std::vector<uint64_t> idxs;
    long double sum_idx = 0.0L;
    long double sum_hash = 0.0L;
    long double sum_log2 = 0.0L;
    long double sum_ratio = 0.0L;
    long double sum_seconds = 0.0L;

    for (const auto& r : all_results) {
        if (r.trunc_bits != t) continue;
        s.repetitions++;
        if (r.found) {
            s.success_count++;
            idxs.push_back(r.collision_index);
            sum_idx += static_cast<long double>(r.collision_index);
            sum_hash += static_cast<long double>(r.hash_evaluations);
            sum_log2 += static_cast<long double>(r.log2_collision_index);
            sum_ratio += r.ratio_to_expected;
            sum_seconds += static_cast<long double>(r.seconds);
        }
    }

    if (idxs.empty()) return s;

    std::sort(idxs.begin(), idxs.end());
    const std::size_t n = idxs.size();
    s.mean_collision_index = sum_idx / static_cast<long double>(n);
    s.mean_hash_evaluations = sum_hash / static_cast<long double>(n);
    s.mean_log2_collision_index = sum_log2 / static_cast<long double>(n);
    s.mean_ratio_to_expected = sum_ratio / static_cast<long double>(n);
    s.mean_seconds = sum_seconds / static_cast<long double>(n);
    s.min_collision_index = idxs.front();
    s.max_collision_index = idxs.back();
    if (n % 2 == 1) s.median_collision_index = static_cast<long double>(idxs[n / 2]);
    else s.median_collision_index = (static_cast<long double>(idxs[n / 2 - 1]) + static_cast<long double>(idxs[n / 2])) / 2.0L;

    if (n > 1) {
        long double variance = 0.0L;
        for (uint64_t x : idxs) {
            const long double diff = static_cast<long double>(x) - s.mean_collision_index;
            variance += diff * diff;
        }
        variance /= static_cast<long double>(n - 1);
        s.stddev_collision_index = std::sqrt(variance);
    }
    return s;
}

static void write_detail_header(std::ofstream& out) {
    out << "algorithm,trunc_bits,repetition,seed_attempt,seed_state_hex,found,"
        << "rho_meet_steps,mu,lambda,collision_index,hash_evaluations,"
        << "expected_collision_steps,log2_collision_index,ratio_to_expected,seconds,"
        << "input_a_hex,input_b_hex,collision_trunc_hex\n";
}

static void write_detail_row(std::ofstream& out, const TrialResult& r) {
    out << "SHA-3-256,"
        << r.trunc_bits << ","
        << r.repetition << ","
        << r.seed_attempt << ","
        << r.seed_state_hex << ","
        << (r.found ? "yes" : "no") << ","
        << r.rho_meet_steps << ","
        << r.mu << ","
        << r.lambda << ","
        << r.collision_index << ","
        << r.hash_evaluations << ","
        << std::fixed << std::setprecision(4) << static_cast<double>(r.expected_steps) << ","
        << std::fixed << std::setprecision(6) << r.log2_collision_index << ","
        << std::fixed << std::setprecision(8) << static_cast<double>(r.ratio_to_expected) << ","
        << std::fixed << std::setprecision(6) << r.seconds << ","
        << r.input_a_hex << ","
        << r.input_b_hex << ","
        << r.collision_trunc_hex << "\n";
}

static void write_summary_csv(const std::string& filepath, const std::vector<SummaryResult>& summaries) {
    std::ofstream out(filepath);
    if (!out) throw std::runtime_error("Failed to open summary CSV output file.");

    out << "algorithm,trunc_bits,repetitions,success_count,expected_collision_steps,"
        << "mean_collision_index,median_collision_index,min_collision_index,max_collision_index,"
        << "stddev_collision_index,mean_hash_evaluations,mean_log2_collision_index,"
        << "mean_ratio_to_expected,mean_seconds\n";

    for (const auto& s : summaries) {
        out << "SHA-3-256,"
            << s.trunc_bits << ","
            << s.repetitions << ","
            << s.success_count << ","
            << std::fixed << std::setprecision(4) << static_cast<double>(s.expected_steps) << ","
            << std::fixed << std::setprecision(4) << static_cast<double>(s.mean_collision_index) << ","
            << std::fixed << std::setprecision(4) << static_cast<double>(s.median_collision_index) << ","
            << s.min_collision_index << ","
            << s.max_collision_index << ","
            << std::fixed << std::setprecision(4) << static_cast<double>(s.stddev_collision_index) << ","
            << std::fixed << std::setprecision(4) << static_cast<double>(s.mean_hash_evaluations) << ","
            << std::fixed << std::setprecision(6) << static_cast<double>(s.mean_log2_collision_index) << ","
            << std::fixed << std::setprecision(8) << static_cast<double>(s.mean_ratio_to_expected) << ","
            << std::fixed << std::setprecision(6) << static_cast<double>(s.mean_seconds) << "\n";
    }
}

static void write_txt_report(const std::string& filepath, const Options& opt, const std::vector<SummaryResult>& summaries) {
    std::ofstream out(filepath);
    if (!out) throw std::runtime_error("Failed to open TXT output file.");

    out << "============================================================\n";
    out << "  MULTICORE TRUNCATED COLLISION EXPERIMENT\n";
    out << "  Algorithm  : SHA-3-256\n";
    out << "============================================================\n\n";
    out << "[Methodology]\n";
    out << "Collision method        : Pollard rho with Floyd tortoise-and-hare\n";
    out << "Projection              : rightmost t bits of 256-bit digest\n";
    out << "Transition              : F_t(x)=rightmost_t(SHA3-256(LE64(salt)||LE64(x)))\n";
    out << "Message length          : 16 bytes\n";
    out << "Maximum phase-1 steps   : " << opt.max_phase_steps << "\n";
    out << "Max seed attempts       : " << opt.max_seed_attempts << "\n";
    out << "Worker budget           : " << opt.workers << "\n";
    out << "Expected baseline       : sqrt(pi/2) * 2^(t/2)\n\n";

    out << "[SHA-3-256 Parameters]\n";
    out << "Permutation             : Keccak-f[1600]\n";
    out << "Rounds                  : 24\n";
    out << "Digest output           : 256 bits\n";
    out << "Capacity                : 512 bits\n";
    out << "Rate                    : 1088 bits\n";
    out << "Seed                    : 0, unseeded mode\n\n";

    out << "[Summary]\n";
    for (const auto& s : summaries) {
        out << "t = " << s.trunc_bits << " bits\n";
        out << "  repetitions              : " << s.repetitions << "\n";
        out << "  success_count            : " << s.success_count << "\n";
        out << "  expected_collision_steps : " << std::fixed << std::setprecision(4) << static_cast<double>(s.expected_steps) << "\n";
        out << "  mean_collision_index     : " << std::fixed << std::setprecision(4) << static_cast<double>(s.mean_collision_index) << "\n";
        out << "  median_collision_index   : " << std::fixed << std::setprecision(4) << static_cast<double>(s.median_collision_index) << "\n";
        out << "  min_collision_index      : " << s.min_collision_index << "\n";
        out << "  max_collision_index      : " << s.max_collision_index << "\n";
        out << "  mean_hash_evaluations    : " << std::fixed << std::setprecision(4) << static_cast<double>(s.mean_hash_evaluations) << "\n";
        out << "  mean_ratio_to_expected   : " << std::fixed << std::setprecision(8) << static_cast<double>(s.mean_ratio_to_expected) << "\n";
        out << "  mean_seconds             : " << std::fixed << std::setprecision(6) << static_cast<double>(s.mean_seconds) << "\n\n";
    }
    out << "=== EXPERIMENT COMPLETED ===\n";
}

static std::vector<TrialResult> run_truncation_level(const Options& opt,
                                                     unsigned t,
                                                     int reps,
                                                     std::ofstream& detail_out,
                                                     std::mutex& detail_mutex,
                                                     std::mutex& cout_mutex) {
    std::vector<TrialResult> results(static_cast<std::size_t>(reps));
    const unsigned controllers_count = std::min<unsigned>(static_cast<unsigned>(reps), std::max(1U, opt.workers));
    std::atomic<int> next_rep(0);
    std::vector<std::thread> controllers;
    controllers.reserve(controllers_count);

    {
        std::lock_guard<std::mutex> lock(cout_mutex);
        std::cout << "[*] trunc_bits=" << t
                  << " repetitions=" << reps
                  << " projection=rightmost"
                  << " parallel_repetitions=" << controllers_count
                  << "\n";
    }

    for (unsigned c = 0; c < controllers_count; ++c) {
        controllers.emplace_back([&, c]() {
            while (true) {
                int rep = next_rep.fetch_add(1, std::memory_order_relaxed);
                if (rep >= reps) break;

                TrialResult r = run_one_trial(t, rep, opt);
                results[static_cast<std::size_t>(rep)] = r;

                {
                    std::lock_guard<std::mutex> lock(detail_mutex);
                    write_detail_row(detail_out, r);
                    detail_out.flush();
                }

                {
                    std::lock_guard<std::mutex> lock(cout_mutex);
                    std::cout << "    rep=" << rep
                              << " found=" << (r.found ? "yes" : "no")
                              << " collision_index=" << r.collision_index
                              << " ratio=" << std::fixed << std::setprecision(6) << static_cast<double>(r.ratio_to_expected)
                              << " evals=" << r.hash_evaluations
                              << " seconds=" << std::fixed << std::setprecision(4) << r.seconds
                              << "\n";
                }
            }
        });
    }

    for (auto& th : controllers) th.join();

    {
        std::lock_guard<std::mutex> lock(cout_mutex);
        std::cout << "\n";
    }

    return results;
}

int main(int argc, char** argv) {
    try {
        Options opt = parse_args(argc, argv);
        const std::string detail_path = opt.out_prefix + "_details.csv";
        const std::string summary_path = opt.out_prefix + "_summary.csv";
        const std::string report_path = opt.out_prefix + "_report.txt";

        std::ofstream detail_out(detail_path);
        if (!detail_out) throw std::runtime_error("Failed to open detailed CSV output file.");
        write_detail_header(detail_out);

        std::mutex detail_mutex;
        std::mutex cout_mutex;
        std::vector<TrialResult> all_results;
        std::vector<SummaryResult> summaries;

        std::cout << "============================================================\n";
        std::cout << "  MULTICORE TRUNCATED COLLISION EXPERIMENT\n";
        std::cout << "  Algorithm  : SHA-3-256\n";
        std::cout << "============================================================\n\n";
        std::cout << "Hardware worker budget : " << opt.workers << "\n";
        std::cout << "Projection             : rightmost t bits\n";
        std::cout << "Message format         : LE64(salt) || LE64(x)\n";
        std::cout << "Max phase-1 steps      : " << opt.max_phase_steps << "\n";
        std::cout << "Max seed attempts      : " << opt.max_seed_attempts << "\n";
        std::cout << "Permutation            : Keccak-f[1600]\n";
        std::cout << "Digest bits            : 256\n\n";

        for (unsigned t : opt.bits) {
            const int reps = repetitions_for_truncation(t);
            std::vector<TrialResult> level_results =
                run_truncation_level(opt, t, reps, detail_out, detail_mutex, cout_mutex);

            all_results.insert(all_results.end(), level_results.begin(), level_results.end());
            summaries.push_back(summarize(t, all_results));
            write_summary_csv(summary_path, summaries);
            write_txt_report(report_path, opt, summaries);
        }

        write_summary_csv(summary_path, summaries);
        write_txt_report(report_path, opt, summaries);

        std::cout << "Done.\n";
        std::cout << "Detailed CSV : " << detail_path << "\n";
        std::cout << "Summary CSV  : " << summary_path << "\n";
        std::cout << "TXT report   : " << report_path << "\n";
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
