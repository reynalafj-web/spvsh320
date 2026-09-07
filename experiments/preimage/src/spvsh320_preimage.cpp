#include <algorithm>
#include <atomic>
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

// ============================================================================
//  MULTICORE UNIFORM TRUNCATED PREIMAGE EXPERIMENT: SpVSH-320
// ============================================================================
//
//  This version is designed for dedicated multicore servers.
//
//  Key changes from the sequential baseline:
//    1. The experiment scenario is aligned with the VSH-256 multicore run.
//    2. Candidate hashing uses a fast 8-byte counter path without repeated
//       vector allocation in the hot loop.
//    3. Repetitions are executed in parallel.
//    4. A single repetition can also use multiple disjoint search lanes when
//       enough hardware threads are available.
//    5. Output is written after each repetition reaches a stable completed state.
//
//  Important methodological note:
//    With lanes_per_rep > 1, candidate search uses disjoint counter chunks across
//    lanes. The reported "trials" value is the total number of candidate hashes
//    actually evaluated by all lanes until the completed state is stable. This is a
//    work-equivalent trial count. It remains comparable to the 2^t geometric
//    expectation because each evaluated candidate is one hash attempt.
//
//  Default experiment:
//    truncation bits : 8,16,24,32,36
//    repetitions     : 50 for 8/16/24/32-bit, 30 for 36-bit
//    max trials      : 2^40 per repetition
//    candidate mode  : counter LE64(i) by default; splitmix64 is available for validation
//    finalization    : 2 extra PermuteSpVSH calls before first output by default
// ============================================================================

// ============================================================================
// SpVSH-320 core implementation, retained from the sequential experiment.
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

static inline uint64_t MulModP1(uint64_t a, uint64_t b) {
    return MulModP(a, b, P1, C1);
}

static inline uint64_t MulModP2(uint64_t a, uint64_t b) {
    return MulModP(a, b, P2, C2);
}

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

static inline uint64_t GET_U64(const uint8_t* ptr) {
    uint64_t val;
    std::memcpy(&val, ptr, 8);
    return val;
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

static inline uint64_t bswap64_portable(uint64_t x) {
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_bswap64(x);
#else
    return ((x & 0x00000000000000FFULL) << 56) |
           ((x & 0x000000000000FF00ULL) << 40) |
           ((x & 0x0000000000FF0000ULL) << 24) |
           ((x & 0x00000000FF000000ULL) << 8)  |
           ((x & 0x000000FF00000000ULL) >> 8)  |
           ((x & 0x0000FF0000000000ULL) >> 24) |
           ((x & 0x00FF000000000000ULL) >> 40) |
           ((x & 0xFF00000000000000ULL) >> 56);
#endif
}

// The sequential baseline serializes the first output lane with memcpy and then
// treats the first 8 digest bytes as a big-endian 64-bit value during truncation.
// On the x86-64 servers used in this experiment, that is equivalent to bswap64(h0).
static inline uint64_t first64_from_output_lane_as_baseline_digest(uint64_t h0) {
    return bswap64_portable(h0);
}

static inline void ApplyFinalizationRounds(PentaState& state, unsigned rounds) {
    for (unsigned r = 0; r < rounds; ++r) {
        PermuteSpVSH(state);
    }
}

// Generic SpVSH path for arbitrary messages. Used for deterministic target-message generation.
static uint64_t hash_spvsh320_first64_message(const uint8_t* in, size_t len, unsigned finalization_rounds) {
    PentaState state = InitSpVSHState();

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

    // Uniform output finalization: apply the same extra permutation depth
    // to every message length before exposing the first squeezed lane.
    ApplyFinalizationRounds(state, finalization_rounds);

    return first64_from_output_lane_as_baseline_digest(state.s[0]);
}

// Fast path for the fixed 8-byte little-endian counter candidate M_i = LE64(i).
static inline uint64_t hash_spvsh320_first64_counter(uint64_t counter, unsigned finalization_rounds) {
    PentaState state = InitSpVSHState();

    // For an 8-byte little-endian candidate, GET_U64(candidate) == counter on x86-64.
    state.s[0] ^= counter;
    PermuteSpVSH(state);

    // Exact 8-byte input leaves an empty final block; the padding block is 0x80 followed by zeros.
    state.s[0] ^= 0x80ULL;
    PermuteSpVSH(state);

    // Uniform output finalization: apply the same extra permutation depth
    // to every message length before exposing the first squeezed lane.
    ApplyFinalizationRounds(state, finalization_rounds);

    return first64_from_output_lane_as_baseline_digest(state.s[0]);
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
        if (v == 0 || v > 63) {
            throw std::runtime_error("truncation bits must be in the range 1..63");
        }
        out.push_back(v);
    }

    if (out.empty()) {
        throw std::runtime_error("empty --bits list");
    }

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
enum class CandidateMode {
    Counter,
    SplitMix64
};

static std::string candidate_mode_to_string(CandidateMode mode) {
    return mode == CandidateMode::SplitMix64 ? "splitmix64_permuted_counter" : "sequential_counter";
}

static CandidateMode parse_candidate_mode(const std::string& s) {
    if (s == "counter" || s == "sequential" || s == "le64") return CandidateMode::Counter;
    if (s == "splitmix64" || s == "splitmix" || s == "mixed" || s == "permuted") return CandidateMode::SplitMix64;
    throw std::runtime_error("invalid candidate mode: " + s + " (use counter or splitmix64)");
}

// Bijective 64-bit mixing function based on the SplitMix64 finalizer.
// Odd multiplications and xor-shift stages make this a deterministic permutation
// over 64-bit words, so the candidate stream is still a one-to-one ordering of
// the 8-byte candidate space rather than a lossy random sample.
static inline uint64_t splitmix64_permute(uint64_t x) {
    x ^= x >> 30;
    x *= UINT64_C(0xBF58476D1CE4E5B9);
    x ^= x >> 27;
    x *= UINT64_C(0x94D049BB133111EB);
    x ^= x >> 31;
    return x;
}

enum class SchedulerMode {
    Fixed,
    Dynamic
};

static std::string scheduler_to_string(SchedulerMode mode) {
    return mode == SchedulerMode::Dynamic ? "dynamic_work_stealing_safe" : "parallel_disjoint_chunked";
}

static SchedulerMode parse_scheduler_mode(const std::string& s) {
    if (s == "fixed" || s == "parallel" || s == "chunked") return SchedulerMode::Fixed;
    if (s == "dynamic" || s == "work-stealing" || s == "work_stealing") return SchedulerMode::Dynamic;
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
    uint64_t trials = 0;              // work-equivalent candidate hashes evaluated
    uint64_t expected_trials = 0;
    uint64_t candidate_counter = 0;   // search index i that produced the match
    uint64_t candidate_word = 0;      // actual 8-byte LE64 candidate word after deterministic mixing
    std::string candidate_mode = "sequential_counter";
    int lanes_per_rep = 1;            // 0 means dynamic worker sharing, not fixed lanes
    std::string search_mode = "parallel_disjoint_chunked";
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
    int lanes_per_rep = 0;             // fixed scheduler only; 0 = auto based on workers/repetitions
    SchedulerMode scheduler = SchedulerMode::Fixed;
    CandidateMode candidate_mode = CandidateMode::Counter;
    uint64_t candidate_seed = UINT64_C(0xD1B54A32D192ED03);
    unsigned finalization_rounds = 2;
    uint64_t max_trials = (uint64_t(1) << 40);
    uint64_t chunk_size = 4096;
    std::string out_prefix = "spvsh320_truncated_preimage_dynamic_safe_finalized";
};

static void print_usage(const char* prog) {
    std::cerr
        << "Usage: " << prog << " [options]\n\n"
        << "Options:\n"
        << "  --bits 8,16,24,32,36       Comma-separated truncation levels\n"
        << "  --workers N                 Total hardware worker budget; 0 = auto\n"
        << "  --scheduler fixed|dynamic   fixed = old bounded lanes, dynamic = safe work stealing\n"
        << "  --candidate-mode MODE       counter or splitmix64; default counter\n"
        << "  --candidate-seed N          64-bit seed/salt for splitmix64 candidate permutation\n"
        << "  --finalization-rounds N     Extra PermuteSpVSH calls before first output; default 2\n"
        << "  --final-rounds N            Alias for --finalization-rounds\n"
        << "  --lanes-per-rep N           Fixed scheduler lanes per repetition; 0 = auto\n"
        << "  --max-trials N              Maximum candidates per repetition\n"
        << "  --max-trials-power P        Set max trials to 2^P\n"
        << "  --chunk-size N              Candidate chunk size per worker, default 4096\n"
        << "  --out-prefix PREFIX         Output filename prefix\n"
        << "  --help                      Show this help\n\n"
        << "Examples:\n"
        << "  " << prog << " --workers 50 --scheduler dynamic --candidate-mode counter --finalization-rounds 2 --bits 8,16,24,32,36 --out-prefix run01_spvsh320_finalized\n"
        << "  " << prog << " --workers 50 --scheduler dynamic --candidate-mode splitmix64 --finalization-rounds 2 --bits 36 --out-prefix validation_mixed_36\n";
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
        } else if (a == "--candidate-mode") {
            opt.candidate_mode = parse_candidate_mode(require_value(a));
        } else if (a == "--candidate-seed") {
            opt.candidate_seed = parse_u64(require_value(a));
        } else if (a == "--finalization-rounds" || a == "--final-rounds") {
            opt.finalization_rounds = static_cast<unsigned>(std::stoul(require_value(a)));
        } else if (a == "--lanes-per-rep") {
            opt.lanes_per_rep = std::stoi(require_value(a));
            if (opt.lanes_per_rep < 0) throw std::runtime_error("--lanes-per-rep must be >= 0");
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

static inline uint64_t repetition_salt(uint64_t seed, int repetition) {
    return seed ^ (UINT64_C(0x9E3779B97F4A7C15) * static_cast<uint64_t>(repetition + 1));
}

static inline uint64_t candidate_word_for_index(uint64_t search_index, int repetition, const Options& opt) {
    if (opt.candidate_mode == CandidateMode::Counter) {
        return search_index;
    }
    const uint64_t salt = repetition_salt(opt.candidate_seed, repetition);
    return splitmix64_permute(search_index ^ salt);
}

// ============================================================================
// Fixed scheduler: bounded repetition-level controllers with optional lanes.
// This preserves the earlier multicore behavior for reproducibility.
// ============================================================================
static TrialResult run_one_trial_fixed(const Options& opt,
                                       unsigned trunc_bits,
                                       int repetition,
                                       uint64_t max_trials,
                                       int lanes_per_rep,
                                       uint64_t chunk_size) {
    const std::string target_message =
        "TRUNCATED_PREIMAGE_TARGET_001|rep=" + std::to_string(repetition);

    const std::vector<uint8_t> target_bytes = string_to_bytes(target_message);
    const uint64_t target_first64 = hash_spvsh320_first64_message(target_bytes.data(), target_bytes.size(), opt.finalization_rounds);
    const uint64_t target_value = leftmost_bits_from_first64(target_first64, trunc_bits);

    TrialResult result;
    result.trunc_bits = trunc_bits;
    result.repetition = repetition;
    result.target_message = target_message;
    result.target_hex = truncated_hex(target_value, trunc_bits);
    result.expected_trials = (uint64_t(1) << trunc_bits);
    result.lanes_per_rep = lanes_per_rep;
    result.search_mode = scheduler_to_string(SchedulerMode::Fixed);
    result.candidate_mode = candidate_mode_to_string(opt.candidate_mode);

    auto start = Clock::now();

    if (lanes_per_rep <= 1) {
        uint64_t attempts = 0;
        for (uint64_t counter = 0; counter < max_trials; ++counter) {
            ++attempts;
            const uint64_t candidate_word = candidate_word_for_index(counter, repetition, opt);
            const uint64_t attempt_first64 = hash_spvsh320_first64_counter(candidate_word, opt.finalization_rounds);
            const uint64_t attempt_value = leftmost_bits_from_first64(attempt_first64, trunc_bits);

            if (attempt_value == target_value) {
                auto end = Clock::now();
                result.found = true;
                result.trials = attempts;
                result.candidate_counter = counter;
                result.candidate_word = candidate_word;
                result.matched_hex = truncated_hex(attempt_value, trunc_bits);
                result.candidate_hex = bytes_to_hex(counter_to_candidate_le(candidate_word));
                result.seconds = std::chrono::duration<double>(end - start).count();
                finalize_trial_metrics(result);
                return result;
            }
        }

        auto end = Clock::now();
        result.found = false;
        result.trials = attempts;
        result.matched_hex = "N/A";
        result.candidate_hex = "N/A";
        result.seconds = std::chrono::duration<double>(end - start).count();
        finalize_trial_metrics(result);
        return result;
    }

    std::atomic<bool> found(false);
    std::atomic<uint64_t> next_counter(0);
    std::mutex winner_mutex;
    std::vector<uint64_t> lane_attempts(static_cast<std::size_t>(lanes_per_rep), 0);

    uint64_t winning_counter = 0;
    uint64_t winning_candidate_word = 0;
    uint64_t winning_value = 0;

    std::vector<std::thread> lanes;
    lanes.reserve(static_cast<std::size_t>(lanes_per_rep));

    for (int lane_id = 0; lane_id < lanes_per_rep; ++lane_id) {
        lanes.emplace_back([&, lane_id]() {
            uint64_t local_attempts = 0;

            while (!found.load(std::memory_order_relaxed)) {
                uint64_t chunk_start = next_counter.fetch_add(chunk_size, std::memory_order_relaxed);
                if (chunk_start >= max_trials) break;

                uint64_t chunk_end = chunk_start + chunk_size;
                if (chunk_end < chunk_start || chunk_end > max_trials) chunk_end = max_trials;

                for (uint64_t counter = chunk_start; counter < chunk_end; ++counter) {
                    if (found.load(std::memory_order_relaxed)) break;

                    ++local_attempts;
                    const uint64_t candidate_word = candidate_word_for_index(counter, repetition, opt);
                    const uint64_t attempt_first64 = hash_spvsh320_first64_counter(candidate_word, opt.finalization_rounds);
                    const uint64_t attempt_value = leftmost_bits_from_first64(attempt_first64, trunc_bits);

                    if (attempt_value == target_value) {
                        bool expected = false;
                        if (found.compare_exchange_strong(expected, true, std::memory_order_relaxed)) {
                            std::lock_guard<std::mutex> lock(winner_mutex);
                            winning_counter = counter;
                            winning_candidate_word = candidate_word;
                            winning_value = attempt_value;
                        }
                        break;
                    }
                }
            }

            lane_attempts[static_cast<std::size_t>(lane_id)] = local_attempts;
        });
    }

    for (auto& th : lanes) th.join();

    auto end = Clock::now();

    uint64_t total_attempts = 0;
    for (uint64_t x : lane_attempts) total_attempts += x;

    result.found = found.load(std::memory_order_relaxed);
    result.trials = total_attempts;
    result.seconds = std::chrono::duration<double>(end - start).count();

    if (result.found) {
        result.candidate_counter = winning_counter;
        result.candidate_word = winning_candidate_word;
        result.matched_hex = truncated_hex(winning_value, trunc_bits);
        result.candidate_hex = bytes_to_hex(counter_to_candidate_le(winning_candidate_word));
    } else {
        result.candidate_counter = 0;
        result.matched_hex = "N/A";
        result.candidate_hex = "N/A";
    }

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
            (static_cast<double>(trials[n / 2 - 1]) +
             static_cast<double>(trials[n / 2])) / 2.0;
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
        << "candidate_length_bytes,candidate_mode,search_mode,lanes_per_rep,found,trials,"
        << "expected_trials,log2_trials,ratio_to_expected,seconds,"
        << "matched_hex,candidate_index,candidate_word,candidate_hex\n";
}

static void write_detail_row(std::ofstream& out, const TrialResult& r) {
    out << "SpVSH-320,"
        << r.trunc_bits << ","
        << r.repetition << ","
        << "\"" << r.target_message << "\","
        << r.target_hex << ","
        << 8 << ","
        << r.candidate_mode << ","
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
        << r.candidate_word << ","
        << r.candidate_hex
        << "\n";
}

static void write_summary_csv(const std::string& filepath,
                              const std::vector<SummaryResult>& summaries) {
    std::ofstream out(filepath);
    if (!out) throw std::runtime_error("Failed to open summary CSV output file.");

    out << "algorithm,trunc_bits,repetitions,success_count,lanes_per_rep,"
        << "expected_trials,mean_trials,median_trials,min_trials,max_trials,"
        << "stddev_trials,mean_log2_trials,mean_ratio_to_expected,mean_seconds\n";

    for (const auto& s : summaries) {
        out << "SpVSH-320,"
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
            << std::fixed << std::setprecision(6) << s.mean_seconds
            << "\n";
    }
}

static void write_txt_report(const std::string& filepath,
                             const Options& opt,
                             const std::vector<SummaryResult>& summaries) {
    std::ofstream out(filepath);
    if (!out) throw std::runtime_error("Failed to open TXT output file.");

    out << "============================================================\n";
    out << "  MULTICORE UNIFORM TRUNCATED PREIMAGE EXPERIMENT\n";
    out << "  Algorithm  : SpVSH-320\n";
    out << "============================================================\n\n";

    out << "[Methodology]\n";
    out << "Target generation       : leftmost t bits of Hash(target_message)\n";
    out << "Candidate generation    : "
        << (opt.candidate_mode == CandidateMode::SplitMix64
            ? "8-byte little-endian encoding of SplitMix64-permuted search index"
            : "8-byte little-endian encoding of counter i") << "\n";
    out << "Candidate mode          : " << candidate_mode_to_string(opt.candidate_mode) << "\n";
    out << "Candidate seed          : 0x" << std::hex << opt.candidate_seed << std::dec << "\n";
    out << "Candidate length        : 8 bytes\n";
    out << "Finalization rounds    : " << opt.finalization_rounds << " extra PermuteSpVSH call(s) before first output\n";
    out << "Scheduler               : " << scheduler_to_string(opt.scheduler) << "\n";
    out << "Search mode             : "
        << (opt.scheduler == SchedulerMode::Dynamic
            ? "dynamic chunk-based work stealing across unfinished repetitions"
            : "parallel disjoint chunked candidate search") << "\n";
    out << "Trial interpretation    : total candidate hashes evaluated until success\n";
    out << "Maximum trials per rep  : " << opt.max_trials << "\n";
    out << "Chunk size              : " << opt.chunk_size << "\n";
    out << "Worker budget           : " << opt.workers << "\n";
    out << "Expected baseline       : 2^t\n\n";

    out << "[SpVSH-320 Parameters]\n";
    out << "Internal state          : 5 x 64-bit lanes = 320 bits\n";
    out << "Digest output           : 256 bits\n";
    out << "Rate lane               : S0, 64 bits\n";
    out << "Moduli                  : P1 = 2^64 - 59, P2 = 2^64 - 83\n";
    out << "Small-prime table       : 64 small odd primes\n\n";

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
    uint64_t winning_candidate_word = 0;
    uint64_t winning_value = 0;
    double seconds = 0.0;
};

static TrialResult dynamic_state_to_result(const DynamicRepState& st, double fallback_seconds, CandidateMode candidate_mode) {
    TrialResult r;
    r.trunc_bits = st.trunc_bits;
    r.repetition = st.repetition;
    r.target_message = st.target_message;
    r.target_hex = st.target_hex;
    r.expected_trials = st.expected_trials;
    r.lanes_per_rep = 0;
    r.search_mode = scheduler_to_string(SchedulerMode::Dynamic);
    r.candidate_mode = candidate_mode_to_string(candidate_mode);
    r.found = st.found.load(std::memory_order_relaxed);
    r.trials = st.attempts.load(std::memory_order_relaxed);
    r.seconds = st.seconds > 0.0 ? st.seconds : fallback_seconds;

    if (r.found) {
        r.candidate_counter = st.winning_counter;
        r.candidate_word = st.winning_candidate_word;
        r.matched_hex = truncated_hex(st.winning_value, st.trunc_bits);
        r.candidate_hex = bytes_to_hex(counter_to_candidate_le(st.winning_candidate_word));
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
        const uint64_t target_first64 = hash_spvsh320_first64_message(target_bytes.data(), target_bytes.size(), opt.finalization_rounds);
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
                  << " scheduler=dynamic_work_stealing_safe"
                  << " workers=" << worker_count
                  << " chunk_size=" << opt.chunk_size
                  << " active_search_lanes=" << worker_count
                  << "\n";
    }

    // Safe finalization rule:
    // A repetition is written only after it has been found or exhausted AND no
    // worker still owns an in-flight chunk for that repetition. This avoids the
    // undercounting that can occur when the winning worker writes the CSV row
    // before other workers have added their local_attempts.
    auto try_finalize = [&](DynamicRepState& st) {
        if (!(st.found.load(std::memory_order_acquire) || st.exhausted.load(std::memory_order_acquire))) {
            return;
        }
        if (st.active_chunks.load(std::memory_order_acquire) != 0) {
            return;
        }

        bool expected_completed = false;
        if (!st.completed.compare_exchange_strong(expected_completed, true, std::memory_order_acq_rel)) {
            return;
        }

        const auto now = Clock::now();
        st.seconds = std::chrono::duration<double>(now - stage_start).count();
        int left = remaining.fetch_sub(1, std::memory_order_acq_rel) - 1;

        TrialResult r = dynamic_state_to_result(st, st.seconds, opt.candidate_mode);
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
            while (remaining.load(std::memory_order_acquire) > 0) {
                bool did_work = false;
                const uint64_t start_idx = cursor.fetch_add(1, std::memory_order_relaxed);

                for (int off = 0; off < reps; ++off) {
                    int idx = static_cast<int>((start_idx + static_cast<uint64_t>(off)) % static_cast<uint64_t>(reps));
                    auto& st = *states[static_cast<std::size_t>(idx)];

                    if (st.completed.load(std::memory_order_acquire) ||
                        st.found.load(std::memory_order_acquire) ||
                        st.exhausted.load(std::memory_order_acquire)) {
                        try_finalize(st);
                        continue;
                    }

                    // Claim an in-flight chunk slot before reading/advancing the
                    // counter. Then re-check the terminal flags to close the race
                    // where another worker finds the target between the first
                    // check and this claim.
                    st.active_chunks.fetch_add(1, std::memory_order_acq_rel);

                    if (st.completed.load(std::memory_order_acquire) ||
                        st.found.load(std::memory_order_acquire) ||
                        st.exhausted.load(std::memory_order_acquire)) {
                        st.active_chunks.fetch_sub(1, std::memory_order_acq_rel);
                        try_finalize(st);
                        continue;
                    }

                    uint64_t chunk_start = st.next_counter.fetch_add(opt.chunk_size, std::memory_order_relaxed);
                    if (chunk_start >= opt.max_trials) {
                        st.exhausted.store(true, std::memory_order_release);
                        st.active_chunks.fetch_sub(1, std::memory_order_acq_rel);
                        try_finalize(st);
                        continue;
                    }

                    did_work = true;
                    uint64_t chunk_end = chunk_start + opt.chunk_size;
                    if (chunk_end < chunk_start || chunk_end > opt.max_trials) chunk_end = opt.max_trials;

                    uint64_t local_attempts = 0;
                    bool matched_here = false;
                    uint64_t local_winning_counter = 0;
                    uint64_t local_winning_candidate_word = 0;
                    uint64_t local_winning_value = 0;

                    for (uint64_t counter = chunk_start; counter < chunk_end; ++counter) {
                        if (st.found.load(std::memory_order_acquire)) break;

                        ++local_attempts;
                        const uint64_t candidate_word = candidate_word_for_index(counter, st.repetition, opt);
                        const uint64_t attempt_first64 = hash_spvsh320_first64_counter(candidate_word, opt.finalization_rounds);
                        const uint64_t attempt_value = leftmost_bits_from_first64(attempt_first64, t);

                        if (attempt_value == st.target_value) {
                            matched_here = true;
                            local_winning_counter = counter;
                            local_winning_candidate_word = candidate_word;
                            local_winning_value = attempt_value;
                            break;
                        }
                    }

                    // Add attempts before publishing the found state, so the
                    // winning worker's own work can never be missing from the
                    // recorded trial count.
                    st.attempts.fetch_add(local_attempts, std::memory_order_acq_rel);

                    if (matched_here) {
                        bool expected_found = false;
                        if (st.found.compare_exchange_strong(expected_found, true, std::memory_order_acq_rel)) {
                            std::lock_guard<std::mutex> lock(st.winner_mutex);
                            st.winning_counter = local_winning_counter;
                            st.winning_candidate_word = local_winning_candidate_word;
                            st.winning_value = local_winning_value;
                        }
                    }

                    st.active_chunks.fetch_sub(1, std::memory_order_acq_rel);
                    try_finalize(st);
                    break;
                }

                if (!did_work) {
                    // Either all repetitions are completed/found/exhausted, or
                    // the only remaining work is currently held by other workers.
                    // Let those workers drain their chunks and finalize safely.
                    for (auto& p : states) {
                        try_finalize(*p);
                    }
                    std::this_thread::yield();
                }
            }
        });
    }

    for (auto& th : workers) th.join();

    const auto stage_end = Clock::now();
    const double stage_seconds = std::chrono::duration<double>(stage_end - stage_start).count();

    // Defensive finalization after all worker threads have joined. At this point
    // no in-flight chunks remain, so every recorded trial count is stable.
    for (auto& p : states) {
        if (!p->completed.load(std::memory_order_acquire)) {
            if (!p->found.load(std::memory_order_acquire) && !p->exhausted.load(std::memory_order_acquire)) {
                p->exhausted.store(true, std::memory_order_release);
            }
            try_finalize(*p);
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
                                                           int lanes_per_rep,
                                                           std::ofstream& detail_out,
                                                           std::mutex& detail_mutex,
                                                           std::mutex& cout_mutex) {
    std::vector<TrialResult> results(static_cast<std::size_t>(reps));

    const unsigned lanes_budget = std::max(1U, opt.workers);
    const unsigned controllers_count = std::min<unsigned>(
        static_cast<unsigned>(reps),
        std::max(1U, lanes_budget / static_cast<unsigned>(lanes_per_rep))
    );
    const unsigned active_search_lanes = controllers_count * static_cast<unsigned>(lanes_per_rep);

    std::atomic<int> next_rep(0);
    std::vector<std::thread> controllers;
    controllers.reserve(controllers_count);

    {
        std::lock_guard<std::mutex> lock(cout_mutex);
        std::cout << "[*] trunc_bits=" << t
                  << " repetitions=" << reps
                  << " scheduler=fixed"
                  << " lanes_per_rep=" << lanes_per_rep
                  << " parallel_repetitions=" << controllers_count
                  << " active_search_lanes=" << active_search_lanes
                  << "\n";
    }

    for (unsigned c = 0; c < controllers_count; ++c) {
        controllers.emplace_back([&, c]() {
            while (true) {
                int rep = next_rep.fetch_add(1, std::memory_order_relaxed);
                if (rep >= reps) break;

                TrialResult r = run_one_trial_fixed(opt, t, rep, opt.max_trials, lanes_per_rep, opt.chunk_size);
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
                              << " counter=" << r.candidate_counter
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
        std::cout << "  Algorithm  : SpVSH-320\n";
        std::cout << "============================================================\n\n";
        std::cout << "Hardware worker budget : " << opt.workers << "\n";
        std::cout << "Scheduler              : " << scheduler_to_string(opt.scheduler) << "\n";
        std::cout << "Candidate mode         : " << candidate_mode_to_string(opt.candidate_mode) << "\n";
        std::cout << "Candidate seed         : 0x" << std::hex << opt.candidate_seed << std::dec << "\n";
        std::cout << "Finalization rounds    : " << opt.finalization_rounds << "\n";
        std::cout << "Max trials per rep     : " << opt.max_trials << "\n";
        std::cout << "Chunk size             : " << opt.chunk_size << "\n";
        std::cout << "SpVSH state bits       : 320\n";
        std::cout << "Digest bits            : 256\n\n";

        for (unsigned t : opt.bits) {
            const int reps = repetitions_for_truncation(t);
            int summary_lanes = 0;
            std::vector<TrialResult> level_results;

            if (opt.scheduler == SchedulerMode::Dynamic) {
                summary_lanes = 0;
                level_results = run_truncation_level_dynamic(opt, t, reps, detail_out, detail_mutex, cout_mutex);
            } else {
                int lanes_per_rep = opt.lanes_per_rep;
                if (lanes_per_rep == 0) {
                    lanes_per_rep = static_cast<int>(std::max(1U, opt.workers / static_cast<unsigned>(reps)));
                }
                if (lanes_per_rep < 1) lanes_per_rep = 1;
                summary_lanes = lanes_per_rep;
                level_results = run_truncation_level_fixed(opt, t, reps, lanes_per_rep,
                                                           detail_out, detail_mutex, cout_mutex);
            }

            all_results.insert(all_results.end(), level_results.begin(), level_results.end());
            summaries.push_back(summarize(t, summary_lanes, all_results));

            // Update partial summary/report after every truncation level.
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
