/*
 * VSH-256 SMHasher3 integration
 *
 * This adapter keeps the VSH-256 construction used in the reduced-target
 * experiments: fixed 128-bit primes p and q, n = p*q, canonical VSH-style
 * block processing, binary length appending, and fixed 256-bit big-endian
 * integer serialization.
 *
 * SMHasher3 passes a seed to every hash function. The VSH-256 construction
 * used in the thesis experiments is unseeded. To satisfy SMHasher3's seed
 * sanity tests without changing the seed=0 mapping, this adapter uses the
 * framework seed only when seed != 0 by domain-prefixing the input with a
 * small SMHasher3-specific seed block. Therefore, seed=0 output matches the
 * VSH-256 experiment implementation exactly, while nonzero seeds are used
 * only for the SMHasher3 test harness.
 */

#include "Platform.h"
#include "Hashlib.h"

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/integer.hpp>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using boost::multiprecision::cpp_int;

// -----------------------------------------------------------------------------
// Fixed 128-bit primes used to construct a 256-bit RSA-style VSH modulus n=p*q.
// These constants match the VSH-256 reduced-target experiment implementation.
// -----------------------------------------------------------------------------
static cpp_int VSH256_fixed_prime_p() {
    return cpp_int("0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF61");
}

static cpp_int VSH256_fixed_prime_q() {
    return cpp_int("0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEED");
}

static cpp_int VSH256_fixed_modulus_n() {
    return VSH256_fixed_prime_p() * VSH256_fixed_prime_q();
}

static unsigned VSH256_bit_length(const cpp_int& x) {
    if (x == 0) return 0;
    return static_cast<unsigned>(boost::multiprecision::msb(x) + 1);
}

static std::vector<uint32_t> VSH256_first_primes(std::size_t count) {
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

static std::vector<int> VSH256_int_to_bits_msb(const cpp_int& x) {
    std::vector<int> bits;
    const unsigned len = VSH256_bit_length(x);
    bits.reserve(len);

    for (int i = static_cast<int>(len) - 1; i >= 0; --i) {
        bits.push_back(static_cast<int>((x >> i) & 1));
    }
    return bits;
}

static std::vector<int> VSH256_encode_length_bits(std::size_t len_value) {
    cpp_int v = len_value;
    if (v == 0) return std::vector<int>{0};
    return VSH256_int_to_bits_msb(v);
}

static std::vector<int> VSH256_bytes_to_bits_msb(const uint8_t* bytes, std::size_t len) {
    std::vector<int> bits;
    bits.reserve(len * 8);

    for (std::size_t bi = 0; bi < len; ++bi) {
        const uint8_t b = bytes[bi];
        for (int shift = 7; shift >= 0; --shift) {
            bits.push_back(static_cast<int>((b >> shift) & 1U));
        }
    }
    return bits;
}

struct VSH256_Context {
    cpp_int n;
    std::vector<uint32_t> primes;
    std::size_t k = 0;
};

static std::size_t VSH256_compute_k_for_modulus(const cpp_int& n,
                                                std::vector<uint32_t>& primes_out) {
    cpp_int prod = 1;
    primes_out.clear();

    const std::vector<uint32_t> primes = VSH256_first_primes(256);
    for (uint32_t p : primes) {
        cpp_int next = prod * p;
        if (next >= n) break;
        prod = next;
        primes_out.push_back(p);
    }
    return primes_out.size();
}

static VSH256_Context VSH256_make_context() {
    VSH256_Context ctx;
    ctx.n = VSH256_fixed_modulus_n();
    ctx.k = VSH256_compute_k_for_modulus(ctx.n, ctx.primes);
    if (ctx.k == 0) {
        throw std::runtime_error("Invalid VSH-256 modulus: computed k = 0");
    }
    return ctx;
}

static const VSH256_Context& VSH256_context() {
    static const VSH256_Context ctx = VSH256_make_context();
    return ctx;
}

// -----------------------------------------------------------------------------
// Canonical VSH-style hash:
//   x_{j+1} = x_j^2 * prod_i p_i^{m_i} mod n
// with zero block padding and binary length appending.
// -----------------------------------------------------------------------------
static cpp_int VSH256_hash_bits(const VSH256_Context& ctx, const std::vector<int>& input_bits) {
    const std::size_t l = input_bits.size();
    const std::size_t L = (l + ctx.k - 1) / ctx.k;
    const std::size_t padded_len = L * ctx.k;

    std::vector<int> m = input_bits;
    m.resize(padded_len, 0);

    std::vector<int> z = VSH256_encode_length_bits(l);
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

static void VSH256_write_32bytes_be(cpp_int x, uint8_t* out) {
    for (int i = 31; i >= 0; --i) {
        out[i] = static_cast<uint8_t>(x & 0xFF);
        x >>= 8;
    }
}

static cpp_int VSH256_hash_message(const uint8_t* data, std::size_t len) {
    const VSH256_Context& ctx = VSH256_context();
    std::vector<int> bits = VSH256_bytes_to_bits_msb(data, len);
    return VSH256_hash_bits(ctx, bits);
}

static void VSH256_digest_unseeded(const uint8_t* data, std::size_t len, uint8_t* out) {
    cpp_int h = VSH256_hash_message(data, len);
    VSH256_write_32bytes_be(h, out);
}

static void VSH256_digest_smhasher_seeded(const uint8_t* data,
                                          std::size_t len,
                                          uint64_t seed,
                                          uint8_t* out) {
    if (seed == 0) {
        VSH256_digest_unseeded(data, len, out);
        return;
    }

    // SMHasher3-only seed adapter. The prefix is applied only for nonzero seed,
    // so seed=0 remains exactly equal to the unseeded VSH-256 construction used
    // in the reduced-target experiments.
    static const uint8_t domain[8] = { 'S', 'M', 'H', '3', 'V', 'S', 'H', 'S' };
    std::vector<uint8_t> msg;
    msg.reserve(16 + len);
    msg.insert(msg.end(), domain, domain + 8);
    for (int i = 0; i < 8; ++i) {
        msg.push_back(static_cast<uint8_t>((seed >> (8 * i)) & 0xFFU));
    }
    msg.insert(msg.end(), data, data + len);

    VSH256_digest_unseeded(msg.data(), msg.size(), out);
}

//------------------------------------------------------------
template <bool bswap>
static void VSH256Hash(const void* in, const size_t len, const seed_t seed, void* out) {
    (void)bswap;
    VSH256_digest_smhasher_seeded(static_cast<const uint8_t*>(in),
                                  len,
                                  static_cast<uint64_t>(seed),
                                  static_cast<uint8_t*>(out));
}

//------------------------------------------------------------
REGISTER_FAMILY(VSH,
   $.src_url    = "local",
   $.src_status = HashFamilyInfo::SRC_UNKNOWN
 );

REGISTER_HASH(VSH_256,
   $.desc            = "VSH-256 fixed-modulus baseline, 256-bit output",
   $.hash_flags      =
         0,
   $.impl_flags      =
         0,
   $.bits            = 256,
   $.verification_LE = 0x0,
   $.verification_BE = 0x0,
   $.hashfn_native   = VSH256Hash<false>,
   $.hashfn_bswap    = VSH256Hash<true>
 );
