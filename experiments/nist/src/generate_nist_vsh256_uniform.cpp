#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/integer.hpp>

#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using boost::multiprecision::cpp_int;

// ============================================================
//  UNIFORM NIST STS DATA GENERATOR: VSH-256
// ============================================================
//
//  Methodology:
//  - Candidate input      : 8-byte little-endian counter LE64(i)
//  - Hash output          : 256-bit digest
//  - Output stream format : ASCII bitstream, characters '0' and '1'
//  - Target size          : 1,000,000,000 bits
//  - NIST STS setting     : 1000 sequences x 1,000,000 bits
//
//  This generator uses the same VSH-256 function family used in the
//  truncated-preimage, truncated-collision, and SMHasher3 integration scripts:
//  - fixed 256-bit RSA-style modulus n = p*q
//  - VSH block processing
//  - zero block padding
//  - binary length appending
//  - 256-bit big-endian digest serialization
// ============================================================

static cpp_int fixed_prime_p() {
    return cpp_int("0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF61");
}

static cpp_int fixed_prime_q() {
    return cpp_int("0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEED");
}

static cpp_int fixed_modulus_n() {
    return fixed_prime_p() * fixed_prime_q();
}

static unsigned bit_length(const cpp_int& x) {
    if (x == 0) return 0;
    return static_cast<unsigned>(boost::multiprecision::msb(x) + 1);
}

static std::vector<uint32_t> first_primes(std::size_t count) {
    std::vector<uint32_t> primes;

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

    if (v == 0) {
        return std::vector<int>(1, 0);
    }

    return int_to_bits_msb(v);
}

static std::vector<int> bytes_to_bits_msb(const uint8_t* data, std::size_t len) {
    std::vector<int> bits;
    bits.reserve(len * 8);

    for (std::size_t i = 0; i < len; ++i) {
        uint8_t b = data[i];

        for (int shift = 7; shift >= 0; --shift) {
            bits.push_back(static_cast<int>((b >> shift) & 1U));
        }
    }

    return bits;
}

static std::vector<uint32_t> compute_primes_for_modulus(const cpp_int& n) {
    cpp_int prod = 1;
    std::vector<uint32_t> primes_out;

    std::vector<uint32_t> primes = first_primes(256);

    for (uint32_t p : primes) {
        cpp_int next = prod * p;

        if (next >= n) break;

        prod = next;
        primes_out.push_back(p);
    }

    return primes_out;
}

struct VSH256Context {
    cpp_int n;
    std::vector<uint32_t> primes;
    std::size_t k;
};

static VSH256Context make_vsh256_context() {
    VSH256Context ctx;
    ctx.n = fixed_modulus_n();
    ctx.primes = compute_primes_for_modulus(ctx.n);
    ctx.k = ctx.primes.size();

    if (ctx.k == 0) {
        throw std::runtime_error("Invalid VSH modulus: computed k = 0.");
    }

    return ctx;
}

static cpp_int vsh_hash_integer(const VSH256Context& ctx,
                                const uint8_t* input,
                                std::size_t len) {
    std::vector<int> input_bits = bytes_to_bits_msb(input, len);

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

static void cpp_int_to_32bytes_be(cpp_int x, uint8_t out[32]) {
    for (int i = 31; i >= 0; --i) {
        out[i] = static_cast<uint8_t>(x & 0xFF);
        x >>= 8;
    }
}

static void hash_vsh256(const VSH256Context& ctx,
                        const uint8_t* input,
                        std::size_t len,
                        uint8_t out[32]) {
    cpp_int h = vsh_hash_integer(ctx, input, len);
    cpp_int_to_32bytes_be(h, out);
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

    std::string filename = "vsh256_nist_1Gbit_ascii.txt";

    if (argc >= 2) {
        filename = argv[1];
    }

    std::cout << "============================================================\n";
    std::cout << "  UNIFORM NIST STS DATA GENERATOR\n";
    std::cout << "  Algorithm  : VSH-256\n";
    std::cout << "============================================================\n\n";

    VSH256Context ctx = make_vsh256_context();

    std::cout << "[INFO] Output file        : " << filename << "\n";
    std::cout << "[INFO] Target bits        : " << TARGET_BITS << "\n";
    std::cout << "[INFO] Digest bits        : " << DIGEST_BITS << "\n";
    std::cout << "[INFO] Iterations         : " << TARGET_ITERATIONS << "\n";
    std::cout << "[INFO] Candidate input    : LE64(counter), 8 bytes\n";
    std::cout << "[INFO] Output format      : ASCII bitstream\n";
    std::cout << "[INFO] NIST STS setting   : 1000 sequences x 1,000,000 bits\n";
    std::cout << "[INFO] VSH modulus bits   : " << bit_length(ctx.n) << "\n";
    std::cout << "[INFO] VSH computed k     : " << ctx.k << "\n\n";

    std::ofstream outfile(filename, std::ios::out | std::ios::binary);

    if (!outfile.is_open()) {
        std::cerr << "[ERROR] Failed to create output file.\n";
        return 1;
    }

    uint8_t payload[8];
    uint8_t digest[32];

    for (uint64_t counter = 0; counter < TARGET_ITERATIONS; ++counter) {
        counter_to_candidate_le(counter, payload);
        hash_vsh256(ctx, payload, 8, digest);
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
