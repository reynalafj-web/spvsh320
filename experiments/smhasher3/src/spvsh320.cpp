/*
 * SpVSH-320-FR2 adapter for SMHasher3
 *
 * This file registers SpVSH-320 as a 256-bit hash in SMHasher3.
 *
 * Important:
 * - The SpVSH-320-FR2 core is kept as in the finalized thesis/paper version:
 *   320-bit state, the same PermuteSpVSH transformation, 0x80 padding,
 *   and two output-finalization rounds before the first output word.
 * - The SMHasher3 seed is used only by this test-harness adapter. When seed = 0,
 *   the output is identical to the unseeded SpVSH-320-FR2 construction used in
 *   the thesis/paper experiments.
 */

#include "Platform.h"
#include "Hashlib.h"

#include <cstdint>
#include <cstring>

//------------------------------------------------------------
// SpVSH-320-FR2 core
//------------------------------------------------------------
struct PentaState {
    uint64_t s[5];
};

static constexpr uint64_t SPVSH_P1 = UINT64_C(0xFFFFFFFFFFFFFFC5); // 2^64 - 59
static constexpr uint64_t SPVSH_C1 = UINT64_C(59);
static constexpr uint64_t SPVSH_P2 = UINT64_C(0xFFFFFFFFFFFFFFAD); // 2^64 - 83
static constexpr uint64_t SPVSH_C2 = UINT64_C(83);
static constexpr unsigned SPVSH_FINALIZATION_ROUNDS = 2;

static const uint16_t SPVSH_SMALL_PRIMES[64] = {
      3,   5,   7,  11,  13,  17,  19,  23,  29,  31,  37,  41,  43,  47,  53,  59,
     61,  67,  71,  73,  79,  83,  89,  97, 101, 103, 107, 109, 113, 127, 131, 137,
    139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227,
    229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313
};

static inline void SpVSH_mul64x64_128(uint64_t a, uint64_t b, uint64_t *res_hi, uint64_t *res_lo) {
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

static inline uint64_t SpVSH_MulModP(uint64_t a, uint64_t b, uint64_t P, uint64_t C) {
    uint64_t hi, lo;
    SpVSH_mul64x64_128(a, b, &hi, &lo);

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

static inline uint64_t SpVSH_MulModP1(uint64_t a, uint64_t b) {
    return SpVSH_MulModP(a, b, SPVSH_P1, SPVSH_C1);
}

static inline uint64_t SpVSH_MulModP2(uint64_t a, uint64_t b) {
    return SpVSH_MulModP(a, b, SPVSH_P2, SPVSH_C2);
}

static inline void SpVSH_Permute(PentaState &state) {
    uint64_t trigger = state.s[0] ^ state.s[2] ^ state.s[4] ^ UINT64_C(0x9E3779B97F4A7C15);

    uint64_t prod_p1 = 1;
    uint64_t prod_p2 = 1;

    while (trigger != 0) {
        int idx = __builtin_ctzll(trigger);
        uint64_t p = SPVSH_SMALL_PRIMES[idx];

        prod_p1 = SpVSH_MulModP1(prod_p1, p);
        prod_p2 = SpVSH_MulModP2(prod_p2, p);

        trigger &= trigger - 1;
    }

    state.s[0] = SpVSH_MulModP1(SpVSH_MulModP1(state.s[0], state.s[0]), prod_p1);
    state.s[1] = SpVSH_MulModP2(SpVSH_MulModP2(state.s[1], state.s[1]), prod_p2);
    state.s[2] = SpVSH_MulModP1(SpVSH_MulModP1(state.s[2], state.s[2]), prod_p1);
    state.s[3] = SpVSH_MulModP2(SpVSH_MulModP2(state.s[3], state.s[3]), prod_p2);
    state.s[4] = SpVSH_MulModP1(SpVSH_MulModP1(state.s[4], state.s[4]), prod_p1);

    uint64_t t0 = state.s[0];
    uint64_t t1 = state.s[1];
    uint64_t t2 = state.s[2];
    uint64_t t3 = state.s[3];
    uint64_t t4 = state.s[4];

    // Same XOR-rotation diffusion layer as the finalized implementation.
    state.s[0] = t0 ^ (t1 >> 17) ^ (t1 << 47);
    state.s[1] = t1 ^ (t2 >> 11) ^ (t2 << 53);
    state.s[2] = t2 ^ (t3 >> 19) ^ (t3 << 45);
    state.s[3] = t3 ^ (t4 >> 13) ^ (t4 << 51);
    state.s[4] = t4 ^ (t0 >> 23) ^ (t0 << 41);
}

static inline PentaState SpVSH_InitState(uint64_t smhasher_seed) {
    PentaState state;
    uint64_t seed64 = UINT64_C(0x0123456789ABCDEF);

    state.s[0] = seed64 + UINT64_C(0xDEADBEEFCAFEBABE);
    state.s[1] = seed64 ^ UINT64_C(0x13198A2E03707344);
    state.s[2] = seed64 + UINT64_C(0x9E3779B97F4A7C15);
    state.s[3] = seed64 ^ UINT64_C(0xBF58476D1CE4E5B9);
    state.s[4] = seed64 + UINT64_C(0x94D049BB133111EB);

    // SMHasher3 adapter seeding.
    // For seed = 0 this has no effect, preserving the paper/thesis SpVSH-320-FR2 output.
    // For seed != 0 this makes the registered SMHasher3 hash seed-sensitive.
    if (smhasher_seed != 0) {
        state.s[1] ^= smhasher_seed;
        state.s[3] ^= smhasher_seed * UINT64_C(0x9E3779B97F4A7C15);
        state.s[4] ^= (smhasher_seed >> 32) | (smhasher_seed << 32);
    }

    for (int i = 0; i < 5; i++) {
        if (state.s[i] == 0) state.s[i] = 3;
    }

    SpVSH_Permute(state);
    return state;
}

static inline void SpVSH_ApplyFinalizationRounds(PentaState& state) {
    for (unsigned r = 0; r < SPVSH_FINALIZATION_ROUNDS; ++r) {
        SpVSH_Permute(state);
    }
}

template <bool bswap>
static void SpVSH320Hash(const void *in, const size_t len, const seed_t seed, void *out) {
    const uint8_t *ptr = static_cast<const uint8_t *>(in);
    size_t remaining = len;
    const uint64_t smhasher_seed = static_cast<uint64_t>(seed);

    PentaState state = SpVSH_InitState(smhasher_seed);

    while (remaining >= 8) {
        uint64_t block = GET_U64<bswap>(ptr, 0);
        state.s[0] ^= block;
        SpVSH_Permute(state);
        ptr += 8;
        remaining -= 8;
    }

    uint8_t buffer[8];
    std::memset(buffer, 0, sizeof(buffer));
    if (remaining > 0) {
        std::memcpy(buffer, ptr, remaining);
    }
    buffer[remaining] = 0x80;

    uint64_t final_block = GET_U64<bswap>(buffer, 0);
    state.s[0] ^= final_block;
    SpVSH_Permute(state);

    SpVSH_ApplyFinalizationRounds(state);

    uint8_t *out8 = static_cast<uint8_t *>(out);
    PUT_U64<bswap>(state.s[0], out8, 0);
    SpVSH_Permute(state);
    PUT_U64<bswap>(state.s[0], out8, 8);
    SpVSH_Permute(state);
    PUT_U64<bswap>(state.s[0], out8, 16);
    SpVSH_Permute(state);
    PUT_U64<bswap>(state.s[0], out8, 24);
}

//------------------------------------------------------------
REGISTER_FAMILY(SpVSH,
   $.src_url    = "local",
   $.src_status = HashFamilyInfo::SRC_UNKNOWN
 );

REGISTER_HASH(SpVSH_320,
   $.desc            = "SpVSH-320-FR2, 256-bit output with SMHasher3 seeded adapter",
   $.hash_flags      =
         0,
   $.impl_flags      =
         0,
   $.bits            = 256,
   $.verification_LE = 0x0,
   $.verification_BE = 0x0,
   $.hashfn_native   = SpVSH320Hash<false>,
   $.hashfn_bswap    = SpVSH320Hash<true>
 );
