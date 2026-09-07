#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/integer.hpp>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
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

using boost::multiprecision::cpp_int;
using Clock = std::chrono::high_resolution_clock;

// ============================================================================
//  MULTICORE TRUNCATED COLLISION EXPERIMENT: VSH-256
// ============================================================================
//
//  Method:
//    Pollard's rho with Floyd's memoryless cycle-finding method.
//
//  Transition function:
//    F_t(x) = rightmost_t(VSH-256(LE64(salt_{t,r}) || LE64(x)))
//
//  Collision extraction:
//    Once a cycle is detected, recover mu and lambda, then verify that
//    F_t(x_{mu-1}) = F_t(x_{mu+lambda-1}) for two distinct predecessor states.
//
//  Default scenario:
//    truncation bits : 16,24,32,40,48,56
//    repetitions     : 50 for 16/24/32-bit, 30 for 40-bit,
//                      10 for 48-bit, 5 for 56-bit
//    phase-1 limit   : 2^32 Floyd iterations
//    projection      : rightmost t bits, i.e., H mod 2^t
//    message format  : LE64(salt_{t,r}) || LE64(x)
//
//  This script preserves the VSH-256 hash implementation style from the
//  previous truncated-preimage experiment. The collision scenario uses a
//  rightmost/modular projection and a 16-byte domain-separated rho message to
//  avoid high-bit fixed-width serialization artifacts.
// ============================================================================

// ----------------------------------------------------------------------------
// Fixed 128-bit primes used to construct a 256-bit RSA-style modulus n = p*q.
// ----------------------------------------------------------------------------
static cpp_int fixed_prime_p() {
    return cpp_int("0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF61");
}

static cpp_int fixed_prime_q() {
    return cpp_int("0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEED");
}

static cpp_int fixed_modulus_n() {
    return fixed_prime_p() * fixed_prime_q();
}

// ----------------------------------------------------------------------------
// Utilities
// ----------------------------------------------------------------------------
static unsigned bit_length(const cpp_int& x) {
    if (x == 0) return 0;
    return static_cast<unsigned>(boost::multiprecision::msb(x) + 1);
}

static std::vector<uint32_t> first_primes(std::size_t count) {
    std::vector<uint32_t> primes;
    primes.reserve(count);

    for (uint32_t x = 2; primes.size() < count; ++x) {
        bool ok = true;
        for (uint32_t p : primes) {
            if (uint64_t(p) * uint64_t(p) > x) break;
            if (x % p == 0) {
                ok = false;
                break;
            }
        }
        if (ok) primes.push_back(x);
    }

    return primes;
}

static std::vector<int> int_to_bits_msb(const cpp_int& x) {
    std::vector<int> bits;
    unsigned len = bit_length(x);
    bits.reserve(len);

    for (int i = static_cast<int>(len) - 1; i >= 0; --i) {
        bits.push_back(static_cast<int>((x >> i) & 1));
    }

    return bits;
}

static std::vector<int> encode_length_bits(std::size_t len_value) {
    cpp_int v = len_value;
    if (v == 0) return {0};
    return int_to_bits_msb(v);
}

static std::vector<int> bytes_to_bits_msb(const std::vector<uint8_t>& bytes) {
    std::vector<int> bits;
    bits.reserve(bytes.size() * 8);

    for (uint8_t b : bytes) {
        for (int shift = 7; shift >= 0; --shift) {
            bits.push_back(static_cast<int>((b >> shift) & 1U));
        }
    }

    return bits;
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

static std::string cpp_int_to_hex(const cpp_int& x) {
    std::ostringstream oss;
    oss << std::hex << x;
    return oss.str();
}

static std::string truncated_hex(uint64_t value, unsigned t) {
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
        if (v == 0 || v > 64) {
            throw std::runtime_error("truncation bits must be in the range 1..64");
        }
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

static uint64_t mask_for_bits(unsigned t) {
    if (t >= 64) return UINT64_MAX;
    return (uint64_t(1) << t) - 1ULL;
}

static inline uint64_t splitmix64_permute(uint64_t x) {
    x += UINT64_C(0x9E3779B97F4A7C15);
    x = (x ^ (x >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
    x = (x ^ (x >> 27)) * UINT64_C(0x94D049BB133111EB);
    x = x ^ (x >> 31);
    return x;
}

static uint64_t deterministic_seed_state(unsigned t, int repetition, int seed_attempt) {
    uint64_t x = UINT64_C(0xD1B54A32D192ED03);
    x ^= UINT64_C(0x9E3779B97F4A7C15) * static_cast<uint64_t>(t + 1);
    x ^= UINT64_C(0xBF58476D1CE4E5B9) * static_cast<uint64_t>(repetition + 1);
    x ^= UINT64_C(0x94D049BB133111EB) * static_cast<uint64_t>(seed_attempt + 1);
    return splitmix64_permute(x) & mask_for_bits(t);
}

static uint64_t deterministic_message_salt(unsigned t, int repetition) {
    uint64_t x = UINT64_C(0xC0111510D15EA5E5);
    x ^= (static_cast<uint64_t>(t) << 32);
    x ^= static_cast<uint64_t>(repetition);
    return splitmix64_permute(x);
}

static std::vector<uint8_t> two_le64_to_message(uint64_t a, uint64_t b) {
    std::vector<uint8_t> out(16, 0);
    for (std::size_t i = 0; i < 8; ++i) out[i] = static_cast<uint8_t>((a >> (8 * i)) & 0xFFU);
    for (std::size_t i = 0; i < 8; ++i) out[8 + i] = static_cast<uint8_t>((b >> (8 * i)) & 0xFFU);
    return out;
}

static double expected_collision_steps(unsigned t) {
    constexpr double sqrt_pi_over_2 = 1.2533141373155002512; // sqrt(pi/2)
    return sqrt_pi_over_2 * std::pow(2.0, static_cast<double>(t) / 2.0);
}

// ----------------------------------------------------------------------------
// VSH context with precomputed parameters.
// ----------------------------------------------------------------------------
struct VSHContext {
    cpp_int n;
    std::vector<uint32_t> primes;
    std::size_t k = 0;

    std::size_t candidate_l = 128;
    std::size_t candidate_L = 0;
    std::size_t candidate_padded_len = 0;
    std::vector<int> candidate_len_bits;
    std::size_t candidate_total_blocks = 0;
};

static std::size_t compute_k_for_modulus(const cpp_int& n,
                                         std::vector<uint32_t>& primes_out) {
    cpp_int prod = 1;
    primes_out.clear();

    auto primes = first_primes(256);
    for (uint32_t p : primes) {
        cpp_int next = prod * p;
        if (next >= n) break;
        prod = next;
        primes_out.push_back(p);
    }

    return primes_out.size();
}

static VSHContext make_vsh_context() {
    VSHContext ctx;
    ctx.n = fixed_modulus_n();
    ctx.k = compute_k_for_modulus(ctx.n, ctx.primes);
    if (ctx.k == 0) throw std::runtime_error("Invalid VSH modulus: computed k = 0.");

    ctx.candidate_l = 128;
    ctx.candidate_L = (ctx.candidate_l + ctx.k - 1) / ctx.k;
    ctx.candidate_padded_len = ctx.candidate_L * ctx.k;
    ctx.candidate_len_bits = encode_length_bits(ctx.candidate_l);
    ctx.candidate_total_blocks =
        (ctx.candidate_padded_len + ctx.candidate_len_bits.size() + ctx.k - 1) / ctx.k;
    return ctx;
}

static uint64_t rightmost_bits_from_cpp_int(const cpp_int& x, unsigned t) {
    if (t == 0 || t > 64) throw std::runtime_error("This implementation supports 1 <= t <= 64.");
    static const cpp_int mask64 = (cpp_int(1) << 64) - 1;
    const uint64_t low64 = (x & mask64).convert_to<uint64_t>();
    if (t == 64) return low64;
    return low64 & ((uint64_t(1) << t) - 1ULL);
}

static cpp_int vsh_hash_bits(const VSHContext& ctx, const std::vector<int>& input_bits) {
    const std::size_t l = input_bits.size();
    const std::size_t L = (l + ctx.k - 1) / ctx.k;
    const std::size_t padded_len = L * ctx.k;

    std::vector<int> m = input_bits;
    m.resize(padded_len, 0);

    std::vector<int> z = encode_length_bits(l);
    m.insert(m.end(), z.begin(), z.end());

    const std::size_t total_blocks = (m.size() + ctx.k - 1) / ctx.k;
    m.resize(total_blocks * ctx.k, 0);

    cpp_int x = 1;
    for (std::size_t j = 0; j < total_blocks; ++j) {
        cpp_int acc = 1;
        for (std::size_t i = 0; i < ctx.k; ++i) {
            if (m[j * ctx.k + i]) acc *= ctx.primes[i];
        }
        x = (x * x) % ctx.n;
        x = (x * acc) % ctx.n;
    }
    return x;
}

static uint64_t hash_vsh256_rightmost_message(const VSHContext& ctx,
                                               const std::vector<uint8_t>& message,
                                               unsigned t) {
    std::vector<int> bits = bytes_to_bits_msb(message);
    cpp_int h = vsh_hash_bits(ctx, bits);
    return rightmost_bits_from_cpp_int(h, t);
}

static inline bool le64_word_bit_msb(uint64_t word, std::size_t bit_index_in_word) {
    // bit_index_in_word follows bytes_to_bits_msb(counter_to_candidate_le(word)).
    // Byte 0 is the least significant byte of the word, but bits inside each byte
    // are read from MSB to LSB.
    const std::size_t byte_index = bit_index_in_word / 8;
    const std::size_t bit_in_byte = bit_index_in_word % 8;
    const unsigned shift = static_cast<unsigned>(8 * byte_index + (7 - bit_in_byte));
    return ((word >> shift) & 1ULL) != 0;
}

static inline bool candidate_16byte_bit_msb(uint64_t salt, uint64_t state, std::size_t bit_index) {
    if (bit_index < 64) {
        return le64_word_bit_msb(salt, bit_index);
    }
    return le64_word_bit_msb(state, bit_index - 64);
}

static uint64_t hash_vsh256_rightmost_salted_state(const VSHContext& ctx,
                                                   uint64_t salt,
                                                   uint64_t state,
                                                   unsigned t) {
    cpp_int x = 1;

    for (std::size_t j = 0; j < ctx.candidate_total_blocks; ++j) {
        cpp_int acc = 1;
        for (std::size_t i = 0; i < ctx.k; ++i) {
            const std::size_t bit_index = j * ctx.k + i;
            bool bit = false;

            if (bit_index < ctx.candidate_l) {
                bit = candidate_16byte_bit_msb(salt, state, bit_index);
            } else if (bit_index >= ctx.candidate_padded_len &&
                       bit_index < ctx.candidate_padded_len + ctx.candidate_len_bits.size()) {
                bit = ctx.candidate_len_bits[bit_index - ctx.candidate_padded_len] != 0;
            }

            if (bit) acc *= ctx.primes[i];
        }
        x = (x * x) % ctx.n;
        x = (x * acc) % ctx.n;
    }

    return rightmost_bits_from_cpp_int(x, t);
}

// ----------------------------------------------------------------------------
// Options and results
// ----------------------------------------------------------------------------
static int repetitions_for_truncation(unsigned t) {
    if (t == 16 || t == 24 || t == 32) return 50;
    if (t == 40) return 30;
    if (t == 48) return 10;
    if (t == 56) return 5;
    return 1;
}

struct Options {
    std::vector<unsigned> bits = {16, 24, 32, 40, 48, 56};
    unsigned workers = 0;
    uint64_t max_phase_steps = (uint64_t(1) << 32);
    int max_seed_attempts = 16;
    std::string out_prefix = "vsh256_truncated_collision_rho";
};

static void print_usage(const char* prog) {
    std::cerr
        << "Usage: " << prog << " [options]\n\n"
        << "Options:\n"
        << "  --bits 16,24,32,40,48,56    Comma-separated truncation levels\n"
        << "  --workers N                  Parallel repetition workers; 0 = auto\n"
        << "  --max-phase N                Maximum Floyd phase-1 iterations per seed attempt\n"
        << "  --max-phase-power P          Set max phase-1 iterations to 2^P\n"
        << "  --max-seed-attempts N        Deterministic seed retries when mu=0 or phase limit is hit\n"
        << "  --out-prefix PREFIX          Output filename prefix\n"
        << "  --help                       Show this help\n\n"
        << "Example:\n"
        << "  " << prog << " --workers 50 --bits 16,24,32,40,48,56 --max-phase-power 32 --out-prefix run01_vsh256_collision\n";
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
            if (opt.max_seed_attempts < 1) throw std::runtime_error("--max-seed-attempts must be >= 1");
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

struct CollisionResult {
    unsigned trunc_bits = 0;
    int repetition = 0;
    int seed_attempt = 0;
    uint64_t seed_state = 0;
    uint64_t message_salt = 0;
    bool found = false;
    uint64_t rho_meet_steps = 0;
    uint64_t mu = 0;
    uint64_t lambda = 0;
    uint64_t collision_index = 0;
    uint64_t hash_evaluations = 0;
    double expected_collision_steps = 0.0;
    double log2_collision_index = 0.0;
    double ratio_to_expected = 0.0;
    double seconds = 0.0;
    uint64_t input_a = 0;
    uint64_t input_b = 0;
    uint64_t collision_value = 0;
    std::string failure_reason;
};

struct SummaryResult {
    unsigned trunc_bits = 0;
    int repetitions = 0;
    int success_count = 0;
    double expected_collision_steps = 0.0;
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

class TransitionFunction {
public:
    TransitionFunction(const VSHContext& ctx, unsigned bits, uint64_t message_salt)
        : ctx_(ctx), bits_(bits), message_salt_(message_salt) {}

    uint64_t operator()(uint64_t x) {
        ++evals_;
        return hash_vsh256_rightmost_salted_state(ctx_, message_salt_, x, bits_);
    }

    uint64_t evals() const { return evals_; }

private:
    const VSHContext& ctx_;
    unsigned bits_;
    uint64_t message_salt_;
    uint64_t evals_ = 0;
};

static CollisionResult run_one_collision_trial(const VSHContext& ctx,
                                               const Options& opt,
                                               unsigned t,
                                               int repetition) {
    CollisionResult result;
    result.trunc_bits = t;
    result.repetition = repetition;
    result.expected_collision_steps = expected_collision_steps(t);

    const auto start_all = Clock::now();
    uint64_t total_evals_before_success = 0;
    const uint64_t message_salt = deterministic_message_salt(t, repetition);
    result.message_salt = message_salt;

    for (int seed_attempt = 0; seed_attempt < opt.max_seed_attempts; ++seed_attempt) {
        TransitionFunction F(ctx, t, message_salt);
        const uint64_t x0 = deterministic_seed_state(t, repetition, seed_attempt);
        result.seed_attempt = seed_attempt;
        result.seed_state = x0;

        uint64_t tortoise = F(x0);
        uint64_t hare = F(F(x0));
        uint64_t meet_steps = 1;
        bool met = false;

        while (meet_steps <= opt.max_phase_steps) {
            if (tortoise == hare) {
                met = true;
                break;
            }
            tortoise = F(tortoise);
            hare = F(F(hare));
            ++meet_steps;
        }

        if (!met) {
            total_evals_before_success += F.evals();
            result.failure_reason = "phase_limit";
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

        if (mu == 0) {
            total_evals_before_success += F.evals();
            result.failure_reason = "mu_zero_retry";
            continue;
        }

        uint64_t input_a = x0;
        for (uint64_t i = 0; i < mu - 1; ++i) {
            input_a = F(input_a);
        }

        uint64_t input_b = x0;
        for (uint64_t i = 0; i < mu + lambda - 1; ++i) {
            input_b = F(input_b);
        }

        const uint64_t fa = F(input_a);
        const uint64_t fb = F(input_b);

        const bool verified = (input_a != input_b) && (fa == fb);
        total_evals_before_success += F.evals();

        if (!verified) {
            result.failure_reason = "verification_failed";
            continue;
        }

        const auto end_all = Clock::now();
        result.found = true;
        result.seed_attempt = seed_attempt;
        result.seed_state = x0;
        result.rho_meet_steps = meet_steps;
        result.mu = mu;
        result.lambda = lambda;
        result.collision_index = mu + lambda;
        result.hash_evaluations = total_evals_before_success;
        result.log2_collision_index = std::log2(static_cast<double>(result.collision_index));
        result.ratio_to_expected = static_cast<double>(result.collision_index) / result.expected_collision_steps;
        result.seconds = std::chrono::duration<double>(end_all - start_all).count();
        result.input_a = input_a;
        result.input_b = input_b;
        result.collision_value = fa;
        result.failure_reason = "";
        return result;
    }

    const auto end_all = Clock::now();
    result.found = false;
    result.hash_evaluations = total_evals_before_success;
    result.seconds = std::chrono::duration<double>(end_all - start_all).count();
    result.collision_index = 0;
    result.log2_collision_index = 0.0;
    result.ratio_to_expected = 0.0;
    if (result.failure_reason.empty()) result.failure_reason = "not_found";
    return result;
}

static SummaryResult summarize(unsigned t, const std::vector<CollisionResult>& all_results) {
    SummaryResult s;
    s.trunc_bits = t;
    s.expected_collision_steps = expected_collision_steps(t);

    std::vector<uint64_t> collision_indices;
    double sum_index = 0.0;
    double sum_hash_evals = 0.0;
    double sum_log2 = 0.0;
    double sum_ratio = 0.0;
    double sum_seconds = 0.0;

    for (const auto& r : all_results) {
        if (r.trunc_bits != t) continue;
        s.repetitions++;
        if (r.found) {
            s.success_count++;
            collision_indices.push_back(r.collision_index);
            sum_index += static_cast<double>(r.collision_index);
            sum_hash_evals += static_cast<double>(r.hash_evaluations);
            sum_log2 += r.log2_collision_index;
            sum_ratio += r.ratio_to_expected;
            sum_seconds += r.seconds;
        }
    }

    if (collision_indices.empty()) return s;

    std::sort(collision_indices.begin(), collision_indices.end());
    const std::size_t n = collision_indices.size();

    s.mean_collision_index = sum_index / static_cast<double>(n);
    s.mean_hash_evaluations = sum_hash_evals / static_cast<double>(n);
    s.mean_log2_collision_index = sum_log2 / static_cast<double>(n);
    s.mean_ratio_to_expected = sum_ratio / static_cast<double>(n);
    s.mean_seconds = sum_seconds / static_cast<double>(n);
    s.min_collision_index = collision_indices.front();
    s.max_collision_index = collision_indices.back();

    if (n % 2 == 1) {
        s.median_collision_index = static_cast<double>(collision_indices[n / 2]);
    } else {
        s.median_collision_index =
            (static_cast<double>(collision_indices[n / 2 - 1]) +
             static_cast<double>(collision_indices[n / 2])) / 2.0;
    }

    if (n > 1) {
        double variance = 0.0;
        for (uint64_t x : collision_indices) {
            const double diff = static_cast<double>(x) - s.mean_collision_index;
            variance += diff * diff;
        }
        variance /= static_cast<double>(n - 1);
        s.stddev_collision_index = std::sqrt(variance);
    }

    return s;
}

// ----------------------------------------------------------------------------
// Output writers
// ----------------------------------------------------------------------------
static void write_detail_header(std::ofstream& out) {
    out << "algorithm,trunc_bits,repetition,seed_attempt,seed_state_hex,message_salt_hex,projection,message_format,found,"
        << "rho_meet_steps,mu,lambda,collision_index,hash_evaluations,"
        << "expected_collision_steps,log2_collision_index,ratio_to_expected,seconds,"
        << "input_a_hex,input_b_hex,collision_trunc_hex,failure_reason\n";
}

static void write_detail_row(std::ofstream& out, const CollisionResult& r) {
    out << "VSH-256,"
        << r.trunc_bits << ","
        << r.repetition << ","
        << r.seed_attempt << ","
        << truncated_hex(r.seed_state, r.trunc_bits) << ","
        << "0x" << std::hex << std::setw(16) << std::setfill('0') << r.message_salt << std::dec << ","
        << "rightmost_mod_2t,"
        << "LE64_salt_concat_LE64_state,"
        << (r.found ? "yes" : "no") << ","
        << r.rho_meet_steps << ","
        << r.mu << ","
        << r.lambda << ","
        << r.collision_index << ","
        << r.hash_evaluations << ","
        << std::fixed << std::setprecision(4) << r.expected_collision_steps << ","
        << std::fixed << std::setprecision(6) << r.log2_collision_index << ","
        << std::fixed << std::setprecision(8) << r.ratio_to_expected << ","
        << std::fixed << std::setprecision(6) << r.seconds << ","
        << (r.found ? bytes_to_hex(two_le64_to_message(r.message_salt, r.input_a)) : "N/A") << ","
        << (r.found ? bytes_to_hex(two_le64_to_message(r.message_salt, r.input_b)) : "N/A") << ","
        << (r.found ? truncated_hex(r.collision_value, r.trunc_bits) : "N/A") << ","
        << (r.failure_reason.empty() ? "" : r.failure_reason)
        << "\n";
}

static void write_summary_csv(const std::string& filepath,
                              const std::vector<SummaryResult>& summaries) {
    std::ofstream out(filepath);
    if (!out) throw std::runtime_error("Failed to open summary CSV output file.");

    out << "algorithm,trunc_bits,repetitions,success_count,expected_collision_steps,"
        << "mean_collision_index,median_collision_index,min_collision_index,max_collision_index,"
        << "stddev_collision_index,mean_hash_evaluations,mean_log2_collision_index,"
        << "mean_ratio_to_expected,mean_seconds\n";

    for (const auto& s : summaries) {
        out << "VSH-256,"
            << s.trunc_bits << ","
            << s.repetitions << ","
            << s.success_count << ","
            << std::fixed << std::setprecision(4) << s.expected_collision_steps << ","
            << std::fixed << std::setprecision(4) << s.mean_collision_index << ","
            << std::fixed << std::setprecision(4) << s.median_collision_index << ","
            << s.min_collision_index << ","
            << s.max_collision_index << ","
            << std::fixed << std::setprecision(4) << s.stddev_collision_index << ","
            << std::fixed << std::setprecision(4) << s.mean_hash_evaluations << ","
            << std::fixed << std::setprecision(6) << s.mean_log2_collision_index << ","
            << std::fixed << std::setprecision(8) << s.mean_ratio_to_expected << ","
            << std::fixed << std::setprecision(6) << s.mean_seconds << "\n";
    }
}

static void write_txt_report(const std::string& filepath,
                             const VSHContext& ctx,
                             const Options& opt,
                             const std::vector<SummaryResult>& summaries) {
    std::ofstream out(filepath);
    if (!out) throw std::runtime_error("Failed to open TXT output file.");

    out << "============================================================\n";
    out << "  MULTICORE TRUNCATED COLLISION EXPERIMENT\n";
    out << "  Algorithm  : VSH-256\n";
    out << "============================================================\n\n";

    out << "[Methodology]\n";
    out << "Collision method        : Pollard rho with Floyd tortoise-and-hare cycle finding\n";
    out << "Transition function     : F_t(x) = rightmost_t(Hash(LE64(salt_{t,r}) || LE64(x)))\n";
    out << "Projection rule         : rightmost t bits, equivalent to H mod 2^t\n";
    out << "Message format          : LE64(salt_{t,r}) || LE64(x)\n";
    out << "Candidate length        : 16 bytes\n";
    out << "Expected baseline       : sqrt(pi/2) * 2^(t/2)\n";
    out << "Max phase-1 iterations  : " << opt.max_phase_steps << "\n";
    out << "Max seed attempts       : " << opt.max_seed_attempts << "\n";
    out << "Worker budget           : " << opt.workers << "\n\n";

    out << "[VSH Parameters]\n";
    out << "p                       : 0x" << cpp_int_to_hex(fixed_prime_p()) << "\n";
    out << "q                       : 0x" << cpp_int_to_hex(fixed_prime_q()) << "\n";
    out << "n = p*q                 : 0x" << cpp_int_to_hex(ctx.n) << "\n";
    out << "Modulus bit-length      : " << bit_length(ctx.n) << " bits\n";
    out << "Computed VSH k          : " << ctx.k << "\n\n";

    out << "[Summary]\n";
    for (const auto& s : summaries) {
        out << "t = " << s.trunc_bits << " bits\n";
        out << "  repetitions             : " << s.repetitions << "\n";
        out << "  success_count           : " << s.success_count << "\n";
        out << "  expected_collision_steps: " << std::fixed << std::setprecision(4) << s.expected_collision_steps << "\n";
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

    out << "=== EXPERIMENT COMPLETED ===\n";
}

static std::vector<CollisionResult> run_truncation_level(const VSHContext& ctx,
                                                         const Options& opt,
                                                         unsigned t,
                                                         int reps,
                                                         std::ofstream& detail_out,
                                                         std::mutex& detail_mutex,
                                                         std::mutex& cout_mutex) {
    std::vector<CollisionResult> results(static_cast<std::size_t>(reps));
    const unsigned worker_count = std::min<unsigned>(static_cast<unsigned>(reps), std::max(1U, opt.workers));
    std::atomic<int> next_rep(0);
    std::vector<std::thread> workers;
    workers.reserve(worker_count);

    {
        std::lock_guard<std::mutex> lock(cout_mutex);
        std::cout << "[*] trunc_bits=" << t
                  << " repetitions=" << reps
                  << " workers=" << worker_count
                  << " expected_steps=" << std::fixed << std::setprecision(2)
                  << expected_collision_steps(t)
                  << "\n";
    }

    for (unsigned w = 0; w < worker_count; ++w) {
        workers.emplace_back([&, w]() {
            while (true) {
                int rep = next_rep.fetch_add(1, std::memory_order_relaxed);
                if (rep >= reps) break;

                CollisionResult r = run_one_collision_trial(ctx, opt, t, rep);
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
                              << " hash_evals=" << r.hash_evaluations
                              << " seconds=" << std::fixed << std::setprecision(4) << r.seconds
                              << (r.failure_reason.empty() ? "" : (" reason=" + r.failure_reason))
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

int main(int argc, char** argv) {
    try {
        Options opt = parse_args(argc, argv);
        VSHContext ctx = make_vsh_context();

        const std::string detail_path = opt.out_prefix + "_details.csv";
        const std::string summary_path = opt.out_prefix + "_summary.csv";
        const std::string report_path = opt.out_prefix + "_report.txt";

        std::ofstream detail_out(detail_path);
        if (!detail_out) throw std::runtime_error("Failed to open detailed CSV output file.");
        write_detail_header(detail_out);

        std::mutex detail_mutex;
        std::mutex cout_mutex;
        std::vector<CollisionResult> all_results;
        std::vector<SummaryResult> summaries;

        std::cout << "============================================================\n";
        std::cout << "  MULTICORE TRUNCATED COLLISION EXPERIMENT\n";
        std::cout << "  Algorithm  : VSH-256\n";
        std::cout << "============================================================\n\n";
        std::cout << "Hardware worker budget : " << opt.workers << "\n";
        std::cout << "Max phase-1 iterations : " << opt.max_phase_steps << "\n";
        std::cout << "Max seed attempts      : " << opt.max_seed_attempts << "\n";
        std::cout << "VSH modulus bits       : " << bit_length(ctx.n) << "\n";
        std::cout << "VSH k                  : " << ctx.k << "\n";
        std::cout << "Projection             : rightmost bits (H mod 2^t)\n";
        std::cout << "Message format         : LE64(salt_{t,r}) || LE64(x)\n\n";

        for (unsigned t : opt.bits) {
            const int reps = repetitions_for_truncation(t);
            std::vector<CollisionResult> level_results =
                run_truncation_level(ctx, opt, t, reps, detail_out, detail_mutex, cout_mutex);

            all_results.insert(all_results.end(), level_results.begin(), level_results.end());
            summaries.push_back(summarize(t, all_results));

            write_summary_csv(summary_path, summaries);
            write_txt_report(report_path, ctx, opt, summaries);
        }

        write_summary_csv(summary_path, summaries);
        write_txt_report(report_path, ctx, opt, summaries);

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
