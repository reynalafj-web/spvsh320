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
//  MULTICORE UNIFORM TRUNCATED PREIMAGE EXPERIMENT: VSH-256
// ============================================================================
//
//  This version is designed for dedicated multicore servers.
//
//  Key changes from the sequential baseline:
//    1. VSH parameters are precomputed once.
//    2. Candidate hashing avoids repeated vector allocation for 8-byte counters.
//    3. Repetitions are executed in parallel.
//    4. A single repetition can also use multiple search lanes when enough
//       hardware threads are available.
//    5. Output is written incrementally after each repetition completes.
//
//  Important methodological note:
//    With lanes_per_rep > 1, candidate search uses disjoint counter chunks across
//    lanes. The reported "trials" value is the total number of candidate hashes
//    actually evaluated by all lanes until one lane finds a match. This is a
//    work-equivalent trial count. It remains comparable to the 2^t geometric
//    expectation because each evaluated candidate is an independent hash attempt.
//
//  Default experiment:
//    truncation bits : 8,16,24,32,36
//    repetitions     : 50 for 8/16/24/32-bit, 30 for 36-bit
//    max trials      : 2^40 per repetition
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

static std::string cpp_int_to_hex(const cpp_int& x) {
    std::ostringstream oss;
    oss << std::hex << x;
    return oss.str();
}

static std::string truncated_hex(uint64_t value, unsigned t) {
    const unsigned width = (t + 3) / 4;

    std::ostringstream oss;
    oss << "0x"
        << std::hex
        << std::setfill('0')
        << std::setw(width)
        << value;

    return oss.str();
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

// ----------------------------------------------------------------------------
// VSH context with precomputed parameters.
// ----------------------------------------------------------------------------
struct VSHContext {
    cpp_int n;
    std::vector<uint32_t> primes;
    std::size_t k = 0;

    // Precomputed layout for the 8-byte counter candidate.
    std::size_t candidate_l = 64;
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

    if (ctx.k == 0) {
        throw std::runtime_error("Invalid VSH modulus: computed k = 0.");
    }

    ctx.candidate_l = 64;
    ctx.candidate_L = (ctx.candidate_l + ctx.k - 1) / ctx.k;
    ctx.candidate_padded_len = ctx.candidate_L * ctx.k;
    ctx.candidate_len_bits = encode_length_bits(ctx.candidate_l);
    ctx.candidate_total_blocks =
        (ctx.candidate_padded_len + ctx.candidate_len_bits.size() + ctx.k - 1) / ctx.k;

    return ctx;
}

static uint64_t first64_from_cpp_int(const cpp_int& x) {
    static const cpp_int mask64 = (cpp_int(1) << 64) - 1;
    cpp_int v = (x >> 192) & mask64;
    return v.convert_to<uint64_t>();
}

static uint64_t leftmost_bits_from_first64(uint64_t first64, unsigned t) {
    if (t == 0 || t > 64) {
        throw std::runtime_error("This implementation supports 1 <= t <= 64.");
    }
    if (t == 64) return first64;
    return first64 >> (64 - t);
}

// ----------------------------------------------------------------------------
// Generic VSH hash for target-message generation. Used only once per repetition.
// ----------------------------------------------------------------------------
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
            if (m[j * ctx.k + i]) {
                acc *= ctx.primes[i];
            }
        }

        x = (x * x) % ctx.n;
        x = (x * acc) % ctx.n;
    }

    return x;
}

static uint64_t hash_vsh256_first64_message(const VSHContext& ctx,
                                            const std::vector<uint8_t>& message) {
    std::vector<int> bits = bytes_to_bits_msb(message);
    cpp_int h = vsh_hash_bits(ctx, bits);
    return first64_from_cpp_int(h);
}

// ----------------------------------------------------------------------------
// Fast VSH path for the 8-byte little-endian counter candidate.
// This avoids repeated vector construction in the hottest loop.
// ----------------------------------------------------------------------------
static inline bool candidate_counter_bit_msb_within_le_bytes(uint64_t counter,
                                                            std::size_t bit_index) {
    // bit_index follows bytes_to_bits_msb(counter_to_candidate_le(counter)).
    // byte 0 is the least significant byte of counter, but bits inside each byte
    // are read from MSB to LSB.
    const std::size_t byte_index = bit_index / 8;
    const std::size_t bit_in_byte = bit_index % 8;
    const unsigned shift = static_cast<unsigned>(8 * byte_index + (7 - bit_in_byte));
    return ((counter >> shift) & 1ULL) != 0;
}

static uint64_t hash_vsh256_first64_counter(const VSHContext& ctx, uint64_t counter) {
    cpp_int x = 1;

    for (std::size_t j = 0; j < ctx.candidate_total_blocks; ++j) {
        cpp_int acc = 1;

        for (std::size_t i = 0; i < ctx.k; ++i) {
            const std::size_t bit_index = j * ctx.k + i;
            bool bit = false;

            if (bit_index < ctx.candidate_l) {
                bit = candidate_counter_bit_msb_within_le_bytes(counter, bit_index);
            } else if (bit_index >= ctx.candidate_padded_len &&
                       bit_index < ctx.candidate_padded_len + ctx.candidate_len_bits.size()) {
                bit = ctx.candidate_len_bits[bit_index - ctx.candidate_padded_len] != 0;
            }

            if (bit) {
                acc *= ctx.primes[i];
            }
        }

        x = (x * x) % ctx.n;
        x = (x * acc) % ctx.n;
    }

    return first64_from_cpp_int(x);
}

// ----------------------------------------------------------------------------
// Experiment data structures
// ----------------------------------------------------------------------------
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
    uint64_t candidate_counter = 0;   // counter value that produced the match
    int lanes_per_rep = 1;
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
    int lanes_per_rep = 0;             // 0 = auto based on workers/repetitions
    uint64_t max_trials = (uint64_t(1) << 40);
    uint64_t chunk_size = 4096;
    std::string out_prefix = "vsh256_truncated_preimage_multicore";
};

static void print_usage(const char* prog) {
    std::cerr
        << "Usage: " << prog << " [options]\n\n"
        << "Options:\n"
        << "  --bits 8,16,24,32,36       Comma-separated truncation levels\n"
        << "  --workers N                 Total hardware worker budget; 0 = auto\n"
        << "  --lanes-per-rep N           Lanes per repetition; 0 = auto\n"
        << "  --max-trials N              Maximum candidates per repetition\n"
        << "  --max-trials-power P        Set max trials to 2^P\n"
        << "  --chunk-size N              Candidate chunk size per lane, default 4096\n"
        << "  --out-prefix PREFIX         Output filename prefix\n"
        << "  --help                      Show this help\n\n"
        << "Example:\n"
        << "  " << prog << " --workers 64 --bits 8,16,24,32,36 --out-prefix run01_vsh256\n";
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

// ----------------------------------------------------------------------------
// One parallel truncated-preimage trial.
// ----------------------------------------------------------------------------
static TrialResult run_one_trial_parallel(const VSHContext& ctx,
                                          unsigned trunc_bits,
                                          int repetition,
                                          uint64_t max_trials,
                                          int lanes_per_rep,
                                          uint64_t chunk_size) {
    const std::string target_message =
        "TRUNCATED_PREIMAGE_TARGET_001|rep=" + std::to_string(repetition);

    const std::vector<uint8_t> target_bytes = string_to_bytes(target_message);
    const uint64_t target_first64 = hash_vsh256_first64_message(ctx, target_bytes);
    const uint64_t target_value = leftmost_bits_from_first64(target_first64, trunc_bits);

    TrialResult result;
    result.trunc_bits = trunc_bits;
    result.repetition = repetition;
    result.target_message = target_message;
    result.target_hex = truncated_hex(target_value, trunc_bits);
    result.expected_trials = (uint64_t(1) << trunc_bits);
    result.lanes_per_rep = lanes_per_rep;

    std::atomic<bool> found(false);
    std::atomic<uint64_t> next_counter(0);
    std::mutex winner_mutex;
    std::vector<uint64_t> lane_attempts(static_cast<std::size_t>(lanes_per_rep), 0);

    uint64_t winning_counter = 0;
    uint64_t winning_value = 0;

    auto start = Clock::now();

    std::vector<std::thread> lanes;
    lanes.reserve(static_cast<std::size_t>(lanes_per_rep));

    for (int lane_id = 0; lane_id < lanes_per_rep; ++lane_id) {
        lanes.emplace_back([&, lane_id]() {
            uint64_t local_attempts = 0;

            while (!found.load(std::memory_order_relaxed)) {
                uint64_t chunk_start = next_counter.fetch_add(chunk_size, std::memory_order_relaxed);
                if (chunk_start >= max_trials) break;

                uint64_t chunk_end = chunk_start + chunk_size;
                if (chunk_end < chunk_start || chunk_end > max_trials) {
                    chunk_end = max_trials;
                }

                for (uint64_t counter = chunk_start; counter < chunk_end; ++counter) {
                    if (found.load(std::memory_order_relaxed)) break;

                    ++local_attempts;
                    const uint64_t attempt_first64 = hash_vsh256_first64_counter(ctx, counter);
                    const uint64_t attempt_value = leftmost_bits_from_first64(attempt_first64, trunc_bits);

                    if (attempt_value == target_value) {
                        bool expected = false;
                        if (found.compare_exchange_strong(expected, true, std::memory_order_relaxed)) {
                            std::lock_guard<std::mutex> lock(winner_mutex);
                            winning_counter = counter;
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
        result.matched_hex = truncated_hex(winning_value, trunc_bits);
        result.candidate_hex = bytes_to_hex(counter_to_candidate_le(winning_counter));
    } else {
        result.candidate_counter = 0;
        result.matched_hex = "N/A";
        result.candidate_hex = "N/A";
    }

    result.log2_trials = result.trials > 0
        ? std::log2(static_cast<double>(result.trials))
        : 0.0;
    result.ratio_to_expected = static_cast<double>(result.trials) /
                               static_cast<double>(result.expected_trials);

    return result;
}

// ----------------------------------------------------------------------------
// Summary computation over successful trials.
// ----------------------------------------------------------------------------
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

// ----------------------------------------------------------------------------
// Output writers
// ----------------------------------------------------------------------------
static void write_detail_header(std::ofstream& out) {
    out << "algorithm,trunc_bits,repetition,target_message,target_hex,"
        << "candidate_length_bytes,search_mode,lanes_per_rep,found,trials,"
        << "expected_trials,log2_trials,ratio_to_expected,seconds,"
        << "matched_hex,candidate_counter,candidate_hex\n";
}

static void write_detail_row(std::ofstream& out, const TrialResult& r) {
    out << "VSH-256,"
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
        out << "VSH-256,"
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
                             const VSHContext& ctx,
                             const Options& opt,
                             const std::vector<SummaryResult>& summaries) {
    std::ofstream out(filepath);
    if (!out) throw std::runtime_error("Failed to open TXT output file.");

    out << "============================================================\n";
    out << "  MULTICORE UNIFORM TRUNCATED PREIMAGE EXPERIMENT\n";
    out << "  Algorithm  : VSH-256\n";
    out << "============================================================\n\n";

    out << "[Methodology]\n";
    out << "Target generation       : leftmost t bits of Hash(target_message)\n";
    out << "Candidate generation    : 8-byte little-endian encoding of counter i\n";
    out << "Candidate length        : 8 bytes\n";
    out << "Search mode             : parallel disjoint chunked candidate search\n";
    out << "Trial interpretation    : total candidate hashes evaluated until success\n";
    out << "Maximum trials per rep  : " << opt.max_trials << "\n";
    out << "Chunk size              : " << opt.chunk_size << "\n";
    out << "Worker budget           : " << opt.workers << "\n";
    out << "Expected baseline       : 2^t\n\n";

    out << "[VSH Parameters]\n";
    out << "p                       : 0x" << cpp_int_to_hex(fixed_prime_p()) << "\n";
    out << "q                       : 0x" << cpp_int_to_hex(fixed_prime_q()) << "\n";
    out << "n = p*q                 : 0x" << cpp_int_to_hex(ctx.n) << "\n";
    out << "Modulus bit-length      : " << bit_length(ctx.n) << " bits\n";
    out << "Computed VSH k          : " << ctx.k << "\n\n";

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

// ----------------------------------------------------------------------------
// Run one truncation level. Repetitions are scheduled with a bounded controller
// pool, so the program does not oversubscribe the server when repetitions exceed
// the requested worker budget.
// ----------------------------------------------------------------------------
static std::vector<TrialResult> run_truncation_level(const VSHContext& ctx,
                                                     const Options& opt,
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

                TrialResult r = run_one_trial_parallel(ctx, t, rep, opt.max_trials,
                                                       lanes_per_rep, opt.chunk_size);
                results[static_cast<std::size_t>(rep)] = r;

                {
                    std::lock_guard<std::mutex> lock(detail_mutex);
                    write_detail_row(detail_out, r);
                    detail_out.flush();
                }

                {
                    std::lock_guard<std::mutex> lock(cout_mutex);
                    const double hps = r.seconds > 0.0
                        ? static_cast<double>(r.trials) / r.seconds
                        : 0.0;

                    std::cout << "    rep=" << rep
                              << " found=" << (r.found ? "yes" : "no")
                              << " trials=" << r.trials
                              << " ratio=" << std::fixed << std::setprecision(6)
                              << r.ratio_to_expected
                              << " seconds=" << std::fixed << std::setprecision(4)
                              << r.seconds
                              << " h/s=" << std::fixed << std::setprecision(2)
                              << hps
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

// ----------------------------------------------------------------------------
// Main
// ----------------------------------------------------------------------------
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

        std::vector<TrialResult> all_results;
        std::vector<SummaryResult> summaries;

        std::cout << "============================================================\n";
        std::cout << "  MULTICORE UNIFORM TRUNCATED PREIMAGE EXPERIMENT\n";
        std::cout << "  Algorithm  : VSH-256\n";
        std::cout << "============================================================\n\n";
        std::cout << "Hardware worker budget : " << opt.workers << "\n";
        std::cout << "Max trials per rep     : " << opt.max_trials << "\n";
        std::cout << "Chunk size             : " << opt.chunk_size << "\n";
        std::cout << "VSH modulus bits       : " << bit_length(ctx.n) << "\n";
        std::cout << "VSH k                  : " << ctx.k << "\n\n";

        for (unsigned t : opt.bits) {
            const int reps = repetitions_for_truncation(t);
            int lanes_per_rep = opt.lanes_per_rep;

            if (lanes_per_rep == 0) {
                lanes_per_rep = static_cast<int>(std::max(1U, opt.workers / static_cast<unsigned>(reps)));
            }
            if (lanes_per_rep < 1) lanes_per_rep = 1;

            std::vector<TrialResult> level_results =
                run_truncation_level(ctx, opt, t, reps, lanes_per_rep,
                                     detail_out, detail_mutex, cout_mutex);

            all_results.insert(all_results.end(), level_results.begin(), level_results.end());
            summaries.push_back(summarize(t, lanes_per_rep, all_results));

            // Update partial summary/report after every truncation level.
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
