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
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

using Clock = std::chrono::high_resolution_clock;

// ============================================================================
//  MULTICORE TRUNCATED COLLISION EXPERIMENT: SpVSH-320-FR2
// ============================================================================
//
//  This script adapts the finalized SpVSH-320 preimage experiment into a
//  truncated-collision experiment. The SpVSH hash core is intentionally kept
//  unchanged from the finalized preimage script: 320-bit state, the same
//  PermuteSpVSH transformation, LE64 absorption, padding, and the default
//  two output-finalization rounds before the first digest word is exposed.
//
//  Collision method:
//    Pollard's rho with Floyd's tortoise-and-hare cycle finding.
//
//  Revised transition function:
//    F_t(x) = rightmost_t( SpVSH-320-FR2( LE64(salt_{t,r}) || LE64(x) ) )
//
//  Rationale for the revised scenario:
//    This aligns the SpVSH-320 collision experiment with the revised VSH-256
//    scenario that uses a domain-separated 16-byte message and a rightmost
//    projection. It avoids making the comparison depend on high-order fixed
//    serialization bits and keeps the same reduced-target collision objective.
//
//  Default truncation levels:
//    16, 24, 32, 40, 48, 56 bits
//
//  Default repetitions:
//    16/24/32-bit : 50 repetitions
//    40-bit       : 30 repetitions
//    48-bit       : 10 repetitions
//    56-bit       : 5 repetitions
//
//  Expected baseline:
//    sqrt(pi/2) * 2^(t/2)
// ============================================================================

// ============================================================================
// SpVSH-320 core implementation, retained from finalized preimage script.
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

static inline uint64_t first64_from_output_lane_as_baseline_digest(uint64_t h0) {
    return bswap64_portable(h0);
}

static inline void ApplyFinalizationRounds(PentaState& state, unsigned rounds) {
    for (unsigned r = 0; r < rounds; ++r) {
        PermuteSpVSH(state);
    }
}

static inline uint64_t output_lane_as_digest_word(uint64_t lane) {
    return first64_from_output_lane_as_baseline_digest(lane);
}

// Fast path for the revised collision message:
//   M_{t,r,x} = LE64(salt_{t,r}) || LE64(x)
// This preserves the SpVSH-320 hash core. Only the experiment message encoding
// changes from the original LE64(x)-only collision pilot.
static inline uint64_t hash_spvsh320_last64_salted_state(uint64_t message_salt,
                                                         uint64_t x,
                                                         unsigned finalization_rounds) {
    PentaState state = InitSpVSHState();

    // First 64-bit block: deterministic domain/repetition salt.
    state.s[0] ^= message_salt;
    PermuteSpVSH(state);

    // Second 64-bit block: rho state x.
    state.s[0] ^= x;
    PermuteSpVSH(state);

    // Exact 16-byte input leaves an empty final block; the padding block is
    // 0x80 followed by zeros.
    state.s[0] ^= 0x80ULL;
    PermuteSpVSH(state);

    // Uniform output finalization retained from the finalized SpVSH-320 design.
    ApplyFinalizationRounds(state, finalization_rounds);

    // Produce the fourth 64-bit digest word H3 and use it as the rightmost
    // 64-bit portion of the 256-bit digest byte string.
    uint64_t h = state.s[0];
    for (int i = 0; i < 3; ++i) {
        PermuteSpVSH(state);
        h = state.s[0];
    }

    return output_lane_as_digest_word(h);
}

// ============================================================================
// Utility functions
// ============================================================================
static uint64_t mask_for_bits(unsigned t) {
    if (t == 0 || t > 64) throw std::runtime_error("t must be in the range 1..64");
    if (t == 64) return UINT64_MAX;
    return (uint64_t(1) << t) - 1ULL;
}

static uint64_t rightmost_bits_from_last64(uint64_t last64, unsigned t) {
    return last64 & mask_for_bits(t);
}

static uint64_t transition_F(uint64_t x,
                             unsigned t,
                             uint64_t message_salt,
                             unsigned finalization_rounds) {
    const uint64_t last64 = hash_spvsh320_last64_salted_state(message_salt, x, finalization_rounds);
    return rightmost_bits_from_last64(last64, t);
}

static std::string u64_hex(uint64_t x, unsigned min_width = 0) {
    std::ostringstream oss;
    oss << "0x" << std::hex << std::setfill('0');
    if (min_width > 0) oss << std::setw(min_width);
    oss << x;
    return oss.str();
}

static std::string le64_hex(uint64_t x) {
    std::ostringstream oss;
    oss << std::hex << std::setfill('0');
    for (int i = 0; i < 8; ++i) {
        oss << std::setw(2) << static_cast<unsigned>((x >> (8 * i)) & 0xFFU);
    }
    return oss.str();
}

static std::string trunc_hex(uint64_t value, unsigned t) {
    const unsigned width = (t + 3) / 4;
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

static inline uint64_t splitmix64_permute(uint64_t x) {
    x += UINT64_C(0x9E3779B97F4A7C15);
    x = (x ^ (x >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
    x = (x ^ (x >> 27)) * UINT64_C(0x94D049BB133111EB);
    x = x ^ (x >> 31);
    return x;
}

static uint64_t seed_state(unsigned t, int repetition, int attempt, uint64_t seed_base) {
    uint64_t x = seed_base;
    x ^= UINT64_C(0xD1B54A32D192ED03) * static_cast<uint64_t>(t + 1);
    x ^= UINT64_C(0x9E3779B97F4A7C15) * static_cast<uint64_t>(repetition + 1);
    x ^= UINT64_C(0xBF58476D1CE4E5B9) * static_cast<uint64_t>(attempt + 1);
    return splitmix64_permute(x) & mask_for_bits(t);
}

static uint64_t deterministic_message_salt(unsigned t, int repetition) {
    uint64_t x = UINT64_C(0xC0111510D15EA5E5);
    x ^= (static_cast<uint64_t>(t) << 32);
    x ^= static_cast<uint64_t>(repetition);
    return splitmix64_permute(x);
}

static std::string le128_hex(uint64_t salt, uint64_t x) {
    std::ostringstream oss;
    oss << std::hex << std::setfill('0');
    for (int i = 0; i < 8; ++i) {
        oss << std::setw(2) << static_cast<unsigned>((salt >> (8 * i)) & 0xFFU);
    }
    for (int i = 0; i < 8; ++i) {
        oss << std::setw(2) << static_cast<unsigned>((x >> (8 * i)) & 0xFFU);
    }
    return oss.str();
}

static int repetitions_for_truncation(unsigned t) {
    if (t == 16 || t == 24 || t == 32) return 50;
    if (t == 40) return 30;
    if (t == 48) return 10;
    if (t == 56) return 5;
    return 10;
}

static double expected_collision_steps(unsigned t) {
    return std::sqrt(M_PI / 2.0) * std::pow(2.0, static_cast<double>(t) / 2.0);
}

// ============================================================================
// Experiment records
// ============================================================================
struct TrialResult {
    unsigned trunc_bits = 0;
    int repetition = 0;
    int seed_attempt = 0;
    uint64_t seed_state_value = 0;
    uint64_t message_salt = 0;
    bool found = false;
    uint64_t rho_meet_steps = 0;
    uint64_t mu = 0;
    uint64_t lambda = 0;
    uint64_t collision_index = 0;
    uint64_t hash_evaluations = 0;
    double expected_steps = 0.0;
    double log2_collision_index = 0.0;
    double ratio_to_expected = 0.0;
    double seconds = 0.0;
    uint64_t input_a = 0;
    uint64_t input_b = 0;
    uint64_t collision_value = 0;
};

struct SummaryResult {
    unsigned trunc_bits = 0;
    int repetitions = 0;
    int success_count = 0;
    double expected_steps = 0.0;
    double mean_collision_index = 0.0;
    double median_collision_index = 0.0;
    uint64_t min_collision_index = 0;
    uint64_t max_collision_index = 0;
    double stddev_collision_index = 0.0;
    double mean_hash_evaluations = 0.0;
    double mean_log2_collision_index = 0.0;
    double mean_ratio_to_expected = 0.0;
    double mean_seconds = 0.0;
};

struct Options {
    std::vector<unsigned> bits = {16, 24, 32, 40, 48, 56};
    unsigned workers = 0;
    unsigned finalization_rounds = 2;
    uint64_t max_phase_steps = (uint64_t(1) << 32);
    int max_seed_attempts = 64;
    uint64_t seed_base = UINT64_C(0xC0111510F00D5EED);
    std::string out_prefix = "spvsh320_truncated_collision_rho_rightproj";
};

static void print_usage(const char* prog) {
    std::cerr
        << "Usage: " << prog << " [options]\n\n"
        << "Options:\n"
        << "  --bits 16,24,32,40,48,56   Comma-separated truncation levels\n"
        << "  --workers N                  Parallel repetition workers; 0 = auto\n"
        << "  --finalization-rounds N      Extra PermuteSpVSH calls before first output; default 2\n"
        << "  --final-rounds N             Alias for --finalization-rounds\n"
        << "  --max-phase-steps N          Maximum Floyd phase-1 iterations\n"
        << "  --max-phase-power P          Set max phase steps to 2^P; default 2^32\n"
        << "  --max-seed-attempts N        Deterministic seed retries when mu=0; default 64\n"
        << "  --seed-base N                64-bit deterministic seed base\n"
        << "  --out-prefix PREFIX          Output filename prefix\n"
        << "  --help                       Show this help\n\n"
        << "Example:\n"
        << "  " << prog << " --workers 50 --finalization-rounds 2 --bits 16,24,32,40,48,56 --out-prefix spvsh320_collision_rightproj_run01\n";
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
        } else if (a == "--finalization-rounds" || a == "--final-rounds") {
            opt.finalization_rounds = static_cast<unsigned>(std::stoul(require_value(a)));
        } else if (a == "--max-phase-steps") {
            opt.max_phase_steps = parse_u64(require_value(a));
        } else if (a == "--max-phase-power") {
            unsigned p = static_cast<unsigned>(std::stoul(require_value(a)));
            if (p >= 63) throw std::runtime_error("--max-phase-power must be < 63");
            opt.max_phase_steps = (uint64_t(1) << p);
        } else if (a == "--max-seed-attempts") {
            opt.max_seed_attempts = std::stoi(require_value(a));
            if (opt.max_seed_attempts <= 0) throw std::runtime_error("--max-seed-attempts must be positive");
        } else if (a == "--seed-base") {
            opt.seed_base = parse_u64(require_value(a));
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

// ============================================================================
// Pollard-rho with Floyd cycle finding
// ============================================================================
static uint64_t iterate_state(uint64_t x0,
                              uint64_t steps,
                              unsigned t,
                              uint64_t message_salt,
                              unsigned finalization_rounds,
                              uint64_t& evals) {
    uint64_t x = x0;
    for (uint64_t i = 0; i < steps; ++i) {
        x = transition_F(x, t, message_salt, finalization_rounds);
        ++evals;
    }
    return x;
}

static TrialResult run_one_trial(const Options& opt, unsigned t, int repetition) {
    const double expected = expected_collision_steps(t);
    const auto trial_start = Clock::now();

    TrialResult result;
    result.trunc_bits = t;
    result.repetition = repetition;
    result.expected_steps = expected;
    result.message_salt = deterministic_message_salt(t, repetition);

    uint64_t total_evals_across_seed_attempts = 0;

    for (int seed_attempt = 0; seed_attempt < opt.max_seed_attempts; ++seed_attempt) {
        uint64_t local_evals = 0;
        const uint64_t x0 = seed_state(t, repetition, seed_attempt, opt.seed_base);

        auto F = [&](uint64_t x) -> uint64_t {
            ++local_evals;
            return transition_F(x, t, result.message_salt, opt.finalization_rounds);
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
            total_evals_across_seed_attempts += local_evals;
            continue;
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

        // If mu=0, there is no predecessor x_{mu-1}. Use another deterministic seed.
        if (mu == 0) {
            total_evals_across_seed_attempts += local_evals;
            continue;
        }

        uint64_t recovery_evals = 0;
        const uint64_t pred_a = iterate_state(x0, mu - 1, t, result.message_salt, opt.finalization_rounds, recovery_evals);
        const uint64_t pred_b = iterate_state(x0, mu + lambda - 1, t, result.message_salt, opt.finalization_rounds, recovery_evals);
        const uint64_t fa = transition_F(pred_a, t, result.message_salt, opt.finalization_rounds); ++recovery_evals;
        const uint64_t fb = transition_F(pred_b, t, result.message_salt, opt.finalization_rounds); ++recovery_evals;

        total_evals_across_seed_attempts += local_evals + recovery_evals;

        const auto trial_end = Clock::now();
        result.seed_attempt = seed_attempt;
        result.seed_state_value = x0;
        result.rho_meet_steps = meet_steps;
        result.mu = mu;
        result.lambda = lambda;
        result.collision_index = mu + lambda;
        result.hash_evaluations = total_evals_across_seed_attempts;
        result.seconds = std::chrono::duration<double>(trial_end - trial_start).count();

        if (fa == fb && pred_a != pred_b) {
            result.found = true;
            result.input_a = pred_a;
            result.input_b = pred_b;
            result.collision_value = fa;
            result.log2_collision_index = std::log2(static_cast<double>(result.collision_index));
            result.ratio_to_expected = static_cast<double>(result.collision_index) / expected;
            return result;
        }
    }

    const auto trial_end = Clock::now();
    result.found = false;
    result.hash_evaluations = total_evals_across_seed_attempts;
    result.seconds = std::chrono::duration<double>(trial_end - trial_start).count();
    result.log2_collision_index = result.collision_index > 0
        ? std::log2(static_cast<double>(result.collision_index))
        : 0.0;
    result.ratio_to_expected = result.collision_index > 0
        ? static_cast<double>(result.collision_index) / expected
        : 0.0;
    return result;
}

// ============================================================================
// Output writers and summary
// ============================================================================
static void write_detail_header(std::ofstream& out) {
    out << "algorithm,trunc_bits,repetition,seed_attempt,seed_state_hex,message_salt_hex,projection,message_format,found,"
        << "rho_meet_steps,mu,lambda,collision_index,hash_evaluations,"
        << "expected_collision_steps,log2_collision_index,ratio_to_expected,seconds,"
        << "input_a_hex,input_b_hex,collision_trunc_hex\n";
}
static void write_detail_row(std::ofstream& out, const TrialResult& r) {
    out << "SpVSH-320," << r.trunc_bits << "," << r.repetition << ","
        << r.seed_attempt << "," << u64_hex(r.seed_state_value, (r.trunc_bits + 3) / 4) << ","
        << "0x" << std::hex << std::setw(16) << std::setfill('0') << r.message_salt << std::dec << ","
        << "rightmost_mod_2t,"
        << "LE64_salt_concat_LE64_state,"
        << (r.found ? "yes" : "no") << ","
        << r.rho_meet_steps << "," << r.mu << "," << r.lambda << ","
        << r.collision_index << "," << r.hash_evaluations << ","
        << std::fixed << std::setprecision(6) << r.expected_steps << ","
        << std::fixed << std::setprecision(6) << r.log2_collision_index << ","
        << std::fixed << std::setprecision(8) << r.ratio_to_expected << ","
        << std::fixed << std::setprecision(6) << r.seconds << ","
        << (r.found ? le128_hex(r.message_salt, r.input_a) : "N/A") << ","
        << (r.found ? le128_hex(r.message_salt, r.input_b) : "N/A") << ","
        << trunc_hex(r.collision_value, r.trunc_bits) << "\n";
}

static SummaryResult summarize(unsigned t, const std::vector<TrialResult>& all_results) {
    SummaryResult s;
    s.trunc_bits = t;
    s.expected_steps = expected_collision_steps(t);

    std::vector<uint64_t> idx;
    double sum_idx = 0.0;
    double sum_eval = 0.0;
    double sum_log2 = 0.0;
    double sum_ratio = 0.0;
    double sum_seconds = 0.0;

    for (const auto& r : all_results) {
        if (r.trunc_bits != t) continue;
        ++s.repetitions;
        if (r.found) {
            ++s.success_count;
            idx.push_back(r.collision_index);
            sum_idx += static_cast<double>(r.collision_index);
            sum_eval += static_cast<double>(r.hash_evaluations);
            sum_log2 += r.log2_collision_index;
            sum_ratio += r.ratio_to_expected;
            sum_seconds += r.seconds;
        }
    }

    if (idx.empty()) return s;
    std::sort(idx.begin(), idx.end());
    const double n = static_cast<double>(idx.size());
    s.mean_collision_index = sum_idx / n;
    s.mean_hash_evaluations = sum_eval / n;
    s.mean_log2_collision_index = sum_log2 / n;
    s.mean_ratio_to_expected = sum_ratio / n;
    s.mean_seconds = sum_seconds / n;
    s.min_collision_index = idx.front();
    s.max_collision_index = idx.back();

    if (idx.size() % 2 == 1) {
        s.median_collision_index = static_cast<double>(idx[idx.size() / 2]);
    } else {
        s.median_collision_index = (static_cast<double>(idx[idx.size() / 2 - 1]) +
                                    static_cast<double>(idx[idx.size() / 2])) / 2.0;
    }

    if (idx.size() > 1) {
        double var = 0.0;
        for (uint64_t x : idx) {
            const double d = static_cast<double>(x) - s.mean_collision_index;
            var += d * d;
        }
        var /= static_cast<double>(idx.size() - 1);
        s.stddev_collision_index = std::sqrt(var);
    }

    return s;
}

static void write_summary_csv(const std::string& filepath, const std::vector<SummaryResult>& summaries) {
    std::ofstream out(filepath);
    if (!out) throw std::runtime_error("Failed to open summary CSV output file.");

    out << "algorithm,trunc_bits,repetitions,success_count,expected_collision_steps,"
        << "mean_collision_index,median_collision_index,min_collision_index,max_collision_index,"
        << "stddev_collision_index,mean_hash_evaluations,mean_log2_collision_index,"
        << "mean_ratio_to_expected,mean_seconds\n";

    for (const auto& s : summaries) {
        out << "SpVSH-320," << s.trunc_bits << "," << s.repetitions << "," << s.success_count << ","
            << std::fixed << std::setprecision(6) << s.expected_steps << ","
            << std::fixed << std::setprecision(4) << s.mean_collision_index << ","
            << std::fixed << std::setprecision(4) << s.median_collision_index << ","
            << s.min_collision_index << "," << s.max_collision_index << ","
            << std::fixed << std::setprecision(4) << s.stddev_collision_index << ","
            << std::fixed << std::setprecision(4) << s.mean_hash_evaluations << ","
            << std::fixed << std::setprecision(6) << s.mean_log2_collision_index << ","
            << std::fixed << std::setprecision(8) << s.mean_ratio_to_expected << ","
            << std::fixed << std::setprecision(6) << s.mean_seconds << "\n";
    }
}

static void write_txt_report(const std::string& filepath, const Options& opt, const std::vector<SummaryResult>& summaries) {
    std::ofstream out(filepath);
    if (!out) throw std::runtime_error("Failed to open TXT report output file.");

    out << "============================================================\n";
    out << "  TRUNCATED COLLISION EXPERIMENT\n";
    out << "  Algorithm  : SpVSH-320\n";
    out << "  Method     : Pollard rho with Floyd cycle finding\n";
    out << "============================================================\n\n";
    out << "[Methodology]\n";
    out << "Transition function    : F_t(x) = rightmost_t(Hash(LE64(salt_{t,r}) || LE64(x)))\n";
    out << "Projection rule        : rightmost t bits, equivalent to H mod 2^t\n";
    out << "Message format         : LE64(salt_{t,r}) || LE64(x)\n";
    out << "Collision extraction   : predecessors x_(mu-1) and x_(mu+lambda-1)\n";
    out << "Expected baseline      : sqrt(pi/2) * 2^(t/2)\n";
    out << "Maximum phase steps    : " << opt.max_phase_steps << "\n";
    out << "Maximum seed attempts  : " << opt.max_seed_attempts << "\n";
    out << "Worker budget          : " << opt.workers << "\n";
    out << "Finalization rounds    : " << opt.finalization_rounds << "\n\n";

    out << "[Summary]\n";
    for (const auto& s : summaries) {
        out << "t = " << s.trunc_bits << " bits\n";
        out << "  repetitions             : " << s.repetitions << "\n";
        out << "  success_count           : " << s.success_count << "\n";
        out << "  expected_collision_steps: " << std::fixed << std::setprecision(6) << s.expected_steps << "\n";
        out << "  mean_collision_index    : " << std::fixed << std::setprecision(4) << s.mean_collision_index << "\n";
        out << "  median_collision_index  : " << std::fixed << std::setprecision(4) << s.median_collision_index << "\n";
        out << "  min_collision_index     : " << s.min_collision_index << "\n";
        out << "  max_collision_index     : " << s.max_collision_index << "\n";
        out << "  stddev_collision_index  : " << std::fixed << std::setprecision(4) << s.stddev_collision_index << "\n";
        out << "  mean_hash_evaluations   : " << std::fixed << std::setprecision(4) << s.mean_hash_evaluations << "\n";
        out << "  mean_log2_collision_idx : " << std::fixed << std::setprecision(6) << s.mean_log2_collision_index << "\n";
        out << "  mean_ratio_to_expected  : " << std::fixed << std::setprecision(8) << s.mean_ratio_to_expected << "\n";
        out << "  mean_seconds            : " << std::fixed << std::setprecision(6) << s.mean_seconds << "\n\n";
    }
}

static std::vector<TrialResult> run_truncation_level(const Options& opt,
                                                     unsigned t,
                                                     int reps,
                                                     std::ofstream& detail_out,
                                                     std::mutex& detail_mutex,
                                                     std::mutex& cout_mutex) {
    std::vector<TrialResult> results(static_cast<std::size_t>(reps));
    std::atomic<int> next_rep(0);
    const unsigned worker_count = std::min<unsigned>(std::max(1U, opt.workers), static_cast<unsigned>(reps));
    std::vector<std::thread> workers;
    workers.reserve(worker_count);

    {
        std::lock_guard<std::mutex> lock(cout_mutex);
        std::cout << "[*] trunc_bits=" << t
                  << " repetitions=" << reps
                  << " workers=" << worker_count
                  << " expected=" << std::fixed << std::setprecision(2) << expected_collision_steps(t)
                  << "\n";
    }

    for (unsigned wid = 0; wid < worker_count; ++wid) {
        workers.emplace_back([&, wid]() {
            while (true) {
                int rep = next_rep.fetch_add(1, std::memory_order_relaxed);
                if (rep >= reps) break;

                TrialResult r = run_one_trial(opt, t, rep);
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
                              << " seed_attempt=" << r.seed_attempt
                              << " collision_index=" << r.collision_index
                              << " ratio=" << std::fixed << std::setprecision(6) << r.ratio_to_expected
                              << " mu=" << r.mu
                              << " lambda=" << r.lambda
                              << " evals=" << r.hash_evaluations
                              << " seconds=" << std::fixed << std::setprecision(4) << r.seconds
                              << "\n";
                }
            }
        });
    }

    for (auto& th : workers) th.join();
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
        std::cout << "  TRUNCATED COLLISION EXPERIMENT\n";
        std::cout << "  Algorithm  : SpVSH-320\n";
        std::cout << "  Method     : Pollard rho / Floyd cycle finding\n";
        std::cout << "============================================================\n\n";
        std::cout << "Hardware worker budget : " << opt.workers << "\n";
        std::cout << "Finalization rounds     : " << opt.finalization_rounds << "\n";
        std::cout << "Max phase steps        : " << opt.max_phase_steps << "\n";
        std::cout << "Max seed attempts      : " << opt.max_seed_attempts << "\n";
        std::cout << "Seed base              : 0x" << std::hex << opt.seed_base << std::dec << "\n\n";

        for (unsigned t : opt.bits) {
            int reps = repetitions_for_truncation(t);
            std::vector<TrialResult> level_results = run_truncation_level(opt, t, reps, detail_out, detail_mutex, cout_mutex);
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
