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
#include <memory>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

using Clock = std::chrono::high_resolution_clock;
using seed_t = uint64_t;

// ============================================================================
//  MULTICORE UNIFORM TRUNCATED PREIMAGE EXPERIMENT: SHA-3-256
// ============================================================================
//
//  Designed for dedicated multicore servers. This version keeps the same
//  target-generation and candidate-generation scenario as the VSH-256 and
//  SpVSH-320 experiments, while supporting both fixed repetition-level parallelism and dynamic work stealing across unfinished repetitions.
//
//  Default experiment:
//    truncation bits : 8,16,24,32,36
//    repetitions     : 50 for 8/16/24/32-bit, 30 for 36-bit
//    max trials      : 2^40 per repetition
//
//  Methodological note:
//    This implementation can use a dynamic chunk-based work-stealing scheduler. Each repetition keeps its own deterministic target and counter space, while worker threads dynamically acquire disjoint counter chunks from unfinished repetitions. Parallelization affects wall-clock time only; the security-oriented metric remains the total number of candidate hashes evaluated until success.
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

/*
 * Kept for compatibility with the uploaded function structure.
 * The experiment uses seed = 0 to keep SHA-3-256 unseeded.
 */
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
    /*
     * SHA-3 domain suffix 01 plus padding.
     */
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

static std::vector<uint8_t> hash_sha3_256_bytes(const std::vector<uint8_t>& message) {
    std::vector<uint8_t> digest(32, 0);

    // seed = 0 preserves unseeded SHA-3 behavior.
    SHA3_256_core<256, false>(message.data(), message.size(), 0, digest.data());

    return digest;
}


// ============================================================================
// Fast first-64-bit extraction paths used by the multicore experiment.
// ============================================================================
static inline uint64_t first64_from_digest_bytes_be_interpretation(const uint8_t* digest) {
    uint64_t first64 = 0;
    for (int i = 0; i < 8; ++i) {
        first64 = (first64 << 8) | digest[i];
    }
    return first64;
}

static uint64_t hash_sha3_256_first64_message(const uint8_t* in, size_t len) {
    uint8_t digest[32] = {0};
    SHA3_256_core<256, false>(in, len, 0, digest);
    return first64_from_digest_bytes_be_interpretation(digest);
}

static inline uint64_t hash_sha3_256_first64_counter(uint64_t counter) {
    sha3_context ctx;
    sha3_Init(&ctx, 256);

    // The baseline candidate message is LE64(counter). For SHA3_Process<false>,
    // absorbing those 8 bytes is equivalent to XORing the counter into lane s[0].
    ctx.s[0] ^= counter;
    ctx.wordIndex = 1;
    ctx.byteIndex = 0;
    ctx.saved = 0;

    uint8_t digest[32] = {0};
    sha3_Finalize<false>(&ctx, 4, digest);
    return first64_from_digest_bytes_be_interpretation(digest);
}

// ============================================================================
// Experiment utilities
// ============================================================================
static std::vector<uint8_t> string_to_bytes(const std::string& s) {
    return std::vector<uint8_t>(s.begin(), s.end());
}

static std::vector<uint8_t> counter_to_candidate_le(uint64_t counter) {
    std::vector<uint8_t> out(8, 0);
    for (std::size_t i = 0; i < 8; ++i) {
        out[i] = static_cast<uint8_t>((counter >> (8 * i)) & 0xFFU);
    }
    return out;
}

static std::string bytes_to_hex(const std::vector<uint8_t>& bytes) {
    std::ostringstream oss;
    oss << std::hex << std::setfill('0');
    for (uint8_t b : bytes) {
        oss << std::setw(2) << static_cast<unsigned>(b);
    }
    return oss.str();
}

static std::string truncated_hex(uint64_t value, unsigned t) {
    const unsigned width = (t + 3) / 4;
    std::ostringstream oss;
    oss << "0x" << std::hex << std::setfill('0') << std::setw(width) << value;
    return oss.str();
}

static uint64_t leftmost_bits_from_first64(uint64_t first64, unsigned t) {
    if (t == 0 || t > 64) {
        throw std::runtime_error("This implementation supports 1 <= t <= 64.");
    }
    if (t == 64) return first64;
    return first64 >> (64 - t);
}

static std::vector<unsigned> parse_bits_list(const std::string& s) {
    std::vector<unsigned> out;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, ',')) {
        if (item.empty()) continue;
        unsigned v = static_cast<unsigned>(std::stoul(item));
        if (v == 0 || v > 63) throw std::runtime_error("truncation bits must be in the range 1..63");
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


// ============================================================================
// Experiment data structures and scheduler options
// ============================================================================
enum class SchedulerMode {
    Fixed,
    Dynamic
};

static std::string scheduler_to_string(SchedulerMode mode) {
    return mode == SchedulerMode::Dynamic ? "dynamic_work_stealing" : "parallel_repetition_level";
}

static SchedulerMode parse_scheduler_mode(const std::string& s) {
    if (s == "fixed" || s == "parallel" || s == "repetition" || s == "repetition-level") {
        return SchedulerMode::Fixed;
    }
    if (s == "dynamic" || s == "work-stealing" || s == "work_stealing") {
        return SchedulerMode::Dynamic;
    }
    throw std::runtime_error("invalid scheduler: " + s + " (use fixed or dynamic)");
}

struct TrialResult {
    unsigned trunc_bits = 0;
    int repetition = 0;
    std::string target_message;
    std::string target_hex;
    std::string matched_hex;
    std::string candidate_hex;
    bool found = false;
    uint64_t trials = 0;              // candidate hashes evaluated for this repetition
    uint64_t expected_trials = 0;
    uint64_t candidate_counter = 0;   // counter value that produced the match
    int lanes_per_rep = 1;            // 0 means dynamic worker sharing
    std::string search_mode = "parallel_repetition_level";
    double log2_trials = 0.0;
    double ratio_to_expected = 0.0;
    double seconds = 0.0;
};

struct SummaryResult {
    unsigned trunc_bits = 0;
    int repetitions = 0;
    int success_count = 0;
    int lanes_per_rep = 1;
    uint64_t expected_trials = 0;
    double mean_trials = 0.0;
    double median_trials = 0.0;
    uint64_t min_trials = 0;
    uint64_t max_trials = 0;
    double stddev_trials = 0.0;
    double mean_log2_trials = 0.0;
    double mean_ratio_to_expected = 0.0;
    double mean_seconds = 0.0;
};

static int repetitions_for_truncation(unsigned t) {
    if (t == 36) return 30;
    return 50;
}

struct Options {
    std::vector<unsigned> bits = {8, 16, 24, 32, 36};
    unsigned workers = 0;              // 0 = hardware_concurrency()
    int lanes_per_rep = 1;             // fixed scheduler only, kept for CLI compatibility
    SchedulerMode scheduler = SchedulerMode::Fixed;
    uint64_t max_trials = (uint64_t(1) << 40);
    uint64_t chunk_size = 4096;
    std::string out_prefix = "sha3_256_truncated_preimage_dynamic";
};

static void print_usage(const char* prog) {
    std::cerr
        << "Usage: " << prog << " [options]\n\n"
        << "Options:\n"
        << "  --bits 8,16,24,32,36       Comma-separated truncation levels\n"
        << "  --workers N                 Total hardware worker budget; 0 = auto\n"
        << "  --scheduler fixed|dynamic   fixed = repetition-level, dynamic = work stealing\n"
        << "  --lanes-per-rep N           Accepted for compatibility; fixed scheduler uses repetition-level workers\n"
        << "  --max-trials N              Maximum candidates per repetition\n"
        << "  --max-trials-power P        Set max trials to 2^P\n"
        << "  --chunk-size N              Candidate chunk size per worker, default 4096\n"
        << "  --out-prefix PREFIX         Output filename prefix\n"
        << "  --help                      Show this help\n\n"
        << "Examples:\n"
        << "  " << prog << " --workers 50 --scheduler dynamic --bits 8,16,24,32,36 --out-prefix run01_sha3\n"
        << "  " << prog << " --workers 50 --scheduler fixed --bits 32,36 --out-prefix fixed_run\n";
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
        } else if (a == "--scheduler") {
            opt.scheduler = parse_scheduler_mode(require_value(a));
        } else if (a == "--lanes-per-rep") {
            opt.lanes_per_rep = std::stoi(require_value(a));
            if (opt.lanes_per_rep < 0) throw std::runtime_error("--lanes-per-rep must be >= 0");
            if (opt.lanes_per_rep == 0) opt.lanes_per_rep = 1;
        } else if (a == "--max-trials") {
            opt.max_trials = parse_u64(require_value(a));
        } else if (a == "--max-trials-power") {
            unsigned p = static_cast<unsigned>(std::stoul(require_value(a)));
            if (p >= 63) throw std::runtime_error("--max-trials-power must be < 63");
            opt.max_trials = (uint64_t(1) << p);
        } else if (a == "--chunk-size") {
            opt.chunk_size = parse_u64(require_value(a));
            if (opt.chunk_size == 0) throw std::runtime_error("--chunk-size must be > 0");
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

static void finalize_trial_metrics(TrialResult& r) {
    r.log2_trials = r.trials > 0 ? std::log2(static_cast<double>(r.trials)) : 0.0;
    r.ratio_to_expected = static_cast<double>(r.trials) / static_cast<double>(r.expected_trials);
}

// ============================================================================
// Fixed scheduler: repetition-level parallelism for reproducibility.
// ============================================================================
static TrialResult run_one_trial_fixed(unsigned trunc_bits,
                                       int repetition,
                                       uint64_t max_trials) {
    const std::string target_message =
        "TRUNCATED_PREIMAGE_TARGET_001|rep=" + std::to_string(repetition);
    const std::vector<uint8_t> target_bytes = string_to_bytes(target_message);
    const uint64_t target_first64 = hash_sha3_256_first64_message(target_bytes.data(), target_bytes.size());
    const uint64_t target_value = leftmost_bits_from_first64(target_first64, trunc_bits);

    TrialResult result;
    result.trunc_bits = trunc_bits;
    result.repetition = repetition;
    result.target_message = target_message;
    result.target_hex = truncated_hex(target_value, trunc_bits);
    result.expected_trials = (uint64_t(1) << trunc_bits);
    result.lanes_per_rep = 1;
    result.search_mode = scheduler_to_string(SchedulerMode::Fixed);

    auto start = Clock::now();
    for (uint64_t counter = 0; counter < max_trials; ++counter) {
        const uint64_t attempt_first64 = hash_sha3_256_first64_counter(counter);
        const uint64_t attempt_value = leftmost_bits_from_first64(attempt_first64, trunc_bits);
        if (attempt_value == target_value) {
            auto end = Clock::now();
            result.found = true;
            result.trials = counter + 1;
            result.candidate_counter = counter;
            result.matched_hex = truncated_hex(attempt_value, trunc_bits);
            result.candidate_hex = bytes_to_hex(counter_to_candidate_le(counter));
            result.seconds = std::chrono::duration<double>(end - start).count();
            finalize_trial_metrics(result);
            return result;
        }
    }

    auto end = Clock::now();
    result.found = false;
    result.trials = max_trials;
    result.matched_hex = "N/A";
    result.candidate_hex = "N/A";
    result.seconds = std::chrono::duration<double>(end - start).count();
    finalize_trial_metrics(result);
    return result;
}

// ============================================================================
// Summary computation over successful trials.
// ============================================================================
static SummaryResult summarize(unsigned trunc_bits,
                               int lanes_per_rep,
                               const std::vector<TrialResult>& all_results) {
    SummaryResult s;
    s.trunc_bits = trunc_bits;
    s.expected_trials = (uint64_t(1) << trunc_bits);
    s.lanes_per_rep = lanes_per_rep;

    std::vector<uint64_t> trials;
    double sum_trials = 0.0;
    double sum_log2 = 0.0;
    double sum_ratio = 0.0;
    double sum_seconds = 0.0;

    for (const auto& r : all_results) {
        if (r.trunc_bits != trunc_bits) continue;
        s.repetitions++;
        if (r.found) {
            s.success_count++;
            trials.push_back(r.trials);
            sum_trials += static_cast<double>(r.trials);
            sum_log2 += r.log2_trials;
            sum_ratio += r.ratio_to_expected;
            sum_seconds += r.seconds;
        }
    }

    if (trials.empty()) return s;
    std::sort(trials.begin(), trials.end());
    s.mean_trials = sum_trials / static_cast<double>(trials.size());
    s.mean_log2_trials = sum_log2 / static_cast<double>(trials.size());
    s.mean_ratio_to_expected = sum_ratio / static_cast<double>(trials.size());
    s.mean_seconds = sum_seconds / static_cast<double>(trials.size());
    s.min_trials = trials.front();
    s.max_trials = trials.back();

    const std::size_t n = trials.size();
    if (n % 2 == 1) {
        s.median_trials = static_cast<double>(trials[n / 2]);
    } else {
        s.median_trials =
            (static_cast<double>(trials[n / 2 - 1]) + static_cast<double>(trials[n / 2])) / 2.0;
    }

    if (trials.size() > 1) {
        double variance = 0.0;
        for (uint64_t x : trials) {
            const double diff = static_cast<double>(x) - s.mean_trials;
            variance += diff * diff;
        }
        variance /= static_cast<double>(trials.size() - 1);
        s.stddev_trials = std::sqrt(variance);
    }
    return s;
}

// ============================================================================
// Output writers
// ============================================================================
static void write_detail_header(std::ofstream& out) {
    out << "algorithm,trunc_bits,repetition,target_message,target_hex,"
        << "candidate_length_bytes,search_mode,lanes_per_rep,found,trials,"
        << "expected_trials,log2_trials,ratio_to_expected,seconds,"
        << "matched_hex,candidate_counter,candidate_hex\n";
}

static void write_detail_row(std::ofstream& out, const TrialResult& r) {
    out << "SHA-3-256,"
        << r.trunc_bits << ","
        << r.repetition << ","
        << "\"" << r.target_message << "\","
        << r.target_hex << ","
        << 8 << ","
        << r.search_mode << ","
        << r.lanes_per_rep << ","
        << (r.found ? "yes" : "no") << ","
        << r.trials << ","
        << r.expected_trials << ","
        << std::fixed << std::setprecision(6) << r.log2_trials << ","
        << std::fixed << std::setprecision(8) << r.ratio_to_expected << ","
        << std::fixed << std::setprecision(6) << r.seconds << ","
        << r.matched_hex << ","
        << r.candidate_counter << ","
        << r.candidate_hex << "\n";
}

static void write_summary_csv(const std::string& filepath,
                              const std::vector<SummaryResult>& summaries) {
    std::ofstream out(filepath);
    if (!out) throw std::runtime_error("Failed to open summary CSV output file.");
    out << "algorithm,trunc_bits,repetitions,success_count,lanes_per_rep,"
        << "expected_trials,mean_trials,median_trials,min_trials,max_trials,"
        << "stddev_trials,mean_log2_trials,mean_ratio_to_expected,mean_seconds\n";
    for (const auto& s : summaries) {
        out << "SHA-3-256,"
            << s.trunc_bits << ","
            << s.repetitions << ","
            << s.success_count << ","
            << s.lanes_per_rep << ","
            << s.expected_trials << ","
            << std::fixed << std::setprecision(4) << s.mean_trials << ","
            << std::fixed << std::setprecision(4) << s.median_trials << ","
            << s.min_trials << ","
            << s.max_trials << ","
            << std::fixed << std::setprecision(4) << s.stddev_trials << ","
            << std::fixed << std::setprecision(6) << s.mean_log2_trials << ","
            << std::fixed << std::setprecision(8) << s.mean_ratio_to_expected << ","
            << std::fixed << std::setprecision(6) << s.mean_seconds << "\n";
    }
}

static void write_txt_report(const std::string& filepath,
                             const Options& opt,
                             const std::vector<SummaryResult>& summaries) {
    std::ofstream out(filepath);
    if (!out) throw std::runtime_error("Failed to open TXT output file.");
    out << "============================================================\n";
    out << "  MULTICORE UNIFORM TRUNCATED PREIMAGE EXPERIMENT\n";
    out << "  Algorithm  : SHA-3-256\n";
    out << "============================================================\n\n";
    out << "[Methodology]\n";
    out << "Target generation       : leftmost t bits of SHA3-256(target_message)\n";
    out << "Candidate generation    : 8-byte little-endian encoding of counter i\n";
    out << "Candidate length        : 8 bytes\n";
    out << "Scheduler               : " << scheduler_to_string(opt.scheduler) << "\n";
    out << "Search mode             : "
        << (opt.scheduler == SchedulerMode::Dynamic
            ? "dynamic chunk-based work stealing across unfinished repetitions"
            : "repetition-level parallelism") << "\n";
    out << "Trial interpretation    : total candidate hashes evaluated until success\n";
    out << "Maximum trials per rep  : " << opt.max_trials << "\n";
    out << "Chunk size              : " << opt.chunk_size << "\n";
    out << "Worker budget           : " << opt.workers << "\n";
    out << "Expected baseline       : 2^t\n\n";
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
        out << "  repetitions           : " << s.repetitions << "\n";
        out << "  success_count         : " << s.success_count << "\n";
        out << "  lanes_per_rep         : " << s.lanes_per_rep << "\n";
        out << "  expected_trials       : " << s.expected_trials << "\n";
        out << "  mean_trials           : " << std::fixed << std::setprecision(4) << s.mean_trials << "\n";
        out << "  median_trials         : " << std::fixed << std::setprecision(4) << s.median_trials << "\n";
        out << "  min_trials            : " << s.min_trials << "\n";
        out << "  max_trials            : " << s.max_trials << "\n";
        out << "  stddev_trials         : " << std::fixed << std::setprecision(4) << s.stddev_trials << "\n";
        out << "  mean_log2_trials      : " << std::fixed << std::setprecision(6) << s.mean_log2_trials << "\n";
        out << "  mean_ratio_to_expected: " << std::fixed << std::setprecision(8) << s.mean_ratio_to_expected << "\n";
        out << "  mean_seconds          : " << std::fixed << std::setprecision(6) << s.mean_seconds << "\n\n";
    }
    out << "=== EXPERIMENT COMPLETED ===\n";
}

// ============================================================================
// Dynamic scheduler state: all worker threads can help any unfinished repetition.
// Each repetition owns an atomic next_counter, so candidate ranges remain disjoint.
// ============================================================================
struct DynamicRepState {
    int repetition = 0;
    unsigned trunc_bits = 0;
    std::string target_message;
    uint64_t target_value = 0;
    std::string target_hex;
    uint64_t expected_trials = 0;

    std::atomic<uint64_t> next_counter{0};
    std::atomic<uint64_t> attempts{0};
    std::atomic<uint32_t> active_chunks{0};
    std::atomic<bool> found{false};
    std::atomic<bool> exhausted{false};
    std::atomic<bool> completed{false};

    std::mutex winner_mutex;
    uint64_t winning_counter = 0;
    uint64_t winning_value = 0;
    double seconds = 0.0;
};

static TrialResult dynamic_state_to_result(const DynamicRepState& st, double fallback_seconds) {
    TrialResult r;
    r.trunc_bits = st.trunc_bits;
    r.repetition = st.repetition;
    r.target_message = st.target_message;
    r.target_hex = st.target_hex;
    r.expected_trials = st.expected_trials;
    r.lanes_per_rep = 0;
    r.search_mode = scheduler_to_string(SchedulerMode::Dynamic);
    r.found = st.found.load(std::memory_order_relaxed);
    r.trials = st.attempts.load(std::memory_order_relaxed);
    r.seconds = st.seconds > 0.0 ? st.seconds : fallback_seconds;

    if (r.found) {
        r.candidate_counter = st.winning_counter;
        r.matched_hex = truncated_hex(st.winning_value, st.trunc_bits);
        r.candidate_hex = bytes_to_hex(counter_to_candidate_le(st.winning_counter));
    } else {
        r.candidate_counter = 0;
        r.matched_hex = "N/A";
        r.candidate_hex = "N/A";
    }
    finalize_trial_metrics(r);
    return r;
}

static std::vector<TrialResult> run_truncation_level_dynamic(const Options& opt,
                                                             unsigned t,
                                                             int reps,
                                                             std::ofstream& detail_out,
                                                             std::mutex& detail_mutex,
                                                             std::mutex& cout_mutex) {
    std::vector<std::shared_ptr<DynamicRepState>> states;
    states.reserve(static_cast<std::size_t>(reps));

    for (int rep = 0; rep < reps; ++rep) {
        auto st = std::make_shared<DynamicRepState>();
        st->repetition = rep;
        st->trunc_bits = t;
        st->target_message = "TRUNCATED_PREIMAGE_TARGET_001|rep=" + std::to_string(rep);
        const std::vector<uint8_t> target_bytes = string_to_bytes(st->target_message);
        const uint64_t target_first64 = hash_sha3_256_first64_message(target_bytes.data(), target_bytes.size());
        st->target_value = leftmost_bits_from_first64(target_first64, t);
        st->target_hex = truncated_hex(st->target_value, t);
        st->expected_trials = (uint64_t(1) << t);
        states.push_back(st);
    }

    std::vector<TrialResult> results(static_cast<std::size_t>(reps));
    std::atomic<int> remaining(reps);
    std::atomic<uint64_t> cursor(0);
    const unsigned worker_count = std::max(1U, opt.workers);
    const auto stage_start = Clock::now();

    {
        std::lock_guard<std::mutex> lock(cout_mutex);
        std::cout << "[*] trunc_bits=" << t
                  << " repetitions=" << reps
                  << " scheduler=dynamic_work_stealing"
                  << " workers=" << worker_count
                  << " chunk_size=" << opt.chunk_size
                  << " active_search_lanes=" << worker_count
                  << "\n";
    }

    auto mark_completed_and_write = [&](DynamicRepState& st) {
        if ((!st.found.load(std::memory_order_relaxed) && !st.exhausted.load(std::memory_order_relaxed)) ||
            st.active_chunks.load(std::memory_order_relaxed) != 0) {
            return;
        }

        bool expected_completed = false;
        if (!st.completed.compare_exchange_strong(expected_completed, true, std::memory_order_relaxed)) {
            return;
        }

        const auto now = Clock::now();
        st.seconds = std::chrono::duration<double>(now - stage_start).count();
        int left = remaining.fetch_sub(1, std::memory_order_relaxed) - 1;

        TrialResult r = dynamic_state_to_result(st, st.seconds);
        results[static_cast<std::size_t>(st.repetition)] = r;

        {
            std::lock_guard<std::mutex> lock(detail_mutex);
            write_detail_row(detail_out, r);
            detail_out.flush();
        }

        {
            std::lock_guard<std::mutex> lock(cout_mutex);
            const double hps = r.seconds > 0.0 ? static_cast<double>(r.trials) / r.seconds : 0.0;
            std::cout << "    rep=" << r.repetition
                      << " found=" << (r.found ? "yes" : "no")
                      << " trials=" << r.trials
                      << " ratio=" << std::fixed << std::setprecision(6) << r.ratio_to_expected
                      << " seconds=" << std::fixed << std::setprecision(4) << r.seconds
                      << " h/s=" << std::fixed << std::setprecision(2) << hps
                      << " counter=" << r.candidate_counter
                      << " remaining=" << left
                      << "\n";
        }
    };

    std::vector<std::thread> workers;
    workers.reserve(worker_count);

    for (unsigned wid = 0; wid < worker_count; ++wid) {
        workers.emplace_back([&, wid]() {
            while (remaining.load(std::memory_order_relaxed) > 0) {
                bool did_work = false;
                const uint64_t start_idx = cursor.fetch_add(1, std::memory_order_relaxed);

                for (int off = 0; off < reps; ++off) {
                    int idx = static_cast<int>((start_idx + static_cast<uint64_t>(off)) % static_cast<uint64_t>(reps));
                    auto& st = *states[static_cast<std::size_t>(idx)];

                    if (st.completed.load(std::memory_order_relaxed)) {
                        continue;
                    }
                    if (st.found.load(std::memory_order_relaxed) || st.exhausted.load(std::memory_order_relaxed)) {
                        mark_completed_and_write(st);
                        continue;
                    }

                    st.active_chunks.fetch_add(1, std::memory_order_relaxed);
                    if (st.completed.load(std::memory_order_relaxed) ||
                        st.found.load(std::memory_order_relaxed) ||
                        st.exhausted.load(std::memory_order_relaxed)) {
                        st.active_chunks.fetch_sub(1, std::memory_order_relaxed);
                        mark_completed_and_write(st);
                        continue;
                    }

                    uint64_t chunk_start = st.next_counter.fetch_add(opt.chunk_size, std::memory_order_relaxed);
                    if (chunk_start >= opt.max_trials) {
                        st.exhausted.store(true, std::memory_order_relaxed);
                        st.active_chunks.fetch_sub(1, std::memory_order_relaxed);
                        mark_completed_and_write(st);
                        continue;
                    }

                    did_work = true;
                    uint64_t chunk_end = chunk_start + opt.chunk_size;
                    if (chunk_end < chunk_start || chunk_end > opt.max_trials) chunk_end = opt.max_trials;

                    uint64_t local_attempts = 0;
                    bool matched_here = false;
                    uint64_t local_winning_counter = 0;
                    uint64_t local_winning_value = 0;

                    for (uint64_t counter = chunk_start; counter < chunk_end; ++counter) {
                        if (st.found.load(std::memory_order_relaxed)) break;
                        ++local_attempts;
                        const uint64_t attempt_first64 = hash_sha3_256_first64_counter(counter);
                        const uint64_t attempt_value = leftmost_bits_from_first64(attempt_first64, t);

                        if (attempt_value == st.target_value) {
                            matched_here = true;
                            local_winning_counter = counter;
                            local_winning_value = attempt_value;
                            break;
                        }
                    }

                    st.attempts.fetch_add(local_attempts, std::memory_order_relaxed);

                    if (matched_here) {
                        bool expected_found = false;
                        if (st.found.compare_exchange_strong(expected_found, true, std::memory_order_relaxed)) {
                            std::lock_guard<std::mutex> lock(st.winner_mutex);
                            st.winning_counter = local_winning_counter;
                            st.winning_value = local_winning_value;
                        }
                    }

                    st.active_chunks.fetch_sub(1, std::memory_order_relaxed);
                    mark_completed_and_write(st);
                    break;
                }

                if (!did_work) {
                    // Try to finalize states that may have become complete while this worker was scanning.
                    for (auto& p : states) {
                        if (!p->completed.load(std::memory_order_relaxed) &&
                            (p->found.load(std::memory_order_relaxed) || p->exhausted.load(std::memory_order_relaxed))) {
                            mark_completed_and_write(*p);
                        }
                    }
                    std::this_thread::yield();
                }
            }
        });
    }

    for (auto& th : workers) th.join();

    const auto stage_end = Clock::now();
    const double stage_seconds = std::chrono::duration<double>(stage_end - stage_start).count();

    for (auto& p : states) {
        if (!p->completed.load(std::memory_order_relaxed)) {
            p->exhausted.store(true, std::memory_order_relaxed);
            mark_completed_and_write(*p);
        }
    }

    {
        std::lock_guard<std::mutex> lock(cout_mutex);
        std::cout << "    dynamic stage wall_seconds=" << std::fixed << std::setprecision(4)
                  << stage_seconds << "\n\n";
    }

    return results;
}

// ============================================================================
// Fixed scheduler truncation level.
// ============================================================================
static std::vector<TrialResult> run_truncation_level_fixed(const Options& opt,
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
                  << " scheduler=fixed"
                  << " parallel_repetitions=" << controllers_count
                  << " search_mode=repetition-level\n";
    }

    for (unsigned c = 0; c < controllers_count; ++c) {
        controllers.emplace_back([&, c]() {
            while (true) {
                int rep = next_rep.fetch_add(1, std::memory_order_relaxed);
                if (rep >= reps) break;
                TrialResult r = run_one_trial_fixed(t, rep, opt.max_trials);
                results[static_cast<std::size_t>(rep)] = r;
                {
                    std::lock_guard<std::mutex> lock(detail_mutex);
                    write_detail_row(detail_out, r);
                    detail_out.flush();
                }
                {
                    std::lock_guard<std::mutex> lock(cout_mutex);
                    const double hps = r.seconds > 0.0 ? static_cast<double>(r.trials) / r.seconds : 0.0;
                    std::cout << "    rep=" << rep
                              << " found=" << (r.found ? "yes" : "no")
                              << " trials=" << r.trials
                              << " ratio=" << std::fixed << std::setprecision(6) << r.ratio_to_expected
                              << " seconds=" << std::fixed << std::setprecision(4) << r.seconds
                              << " h/s=" << std::fixed << std::setprecision(2) << hps
                              << " counter=" << r.candidate_counter << "\n";
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

// ============================================================================
// Main
// ============================================================================
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
        std::cout << "  MULTICORE UNIFORM TRUNCATED PREIMAGE EXPERIMENT\n";
        std::cout << "  Algorithm  : SHA-3-256\n";
        std::cout << "============================================================\n\n";
        std::cout << "Hardware worker budget : " << opt.workers << "\n";
        std::cout << "Scheduler              : " << scheduler_to_string(opt.scheduler) << "\n";
        std::cout << "Max trials per rep     : " << opt.max_trials << "\n";
        std::cout << "Chunk size             : " << opt.chunk_size << "\n";
        std::cout << "Permutation            : Keccak-f[1600]\n";
        std::cout << "Digest bits            : 256\n";
        std::cout << "Seed                   : 0, unseeded mode\n\n";

        for (unsigned t : opt.bits) {
            const int reps = repetitions_for_truncation(t);
            int summary_lanes = 0;
            std::vector<TrialResult> level_results;

            if (opt.scheduler == SchedulerMode::Dynamic) {
                summary_lanes = 0;
                level_results = run_truncation_level_dynamic(opt, t, reps, detail_out, detail_mutex, cout_mutex);
            } else {
                summary_lanes = 1;
                level_results = run_truncation_level_fixed(opt, t, reps, detail_out, detail_mutex, cout_mutex);
            }

            all_results.insert(all_results.end(), level_results.begin(), level_results.end());
            summaries.push_back(summarize(t, summary_lanes, all_results));
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
