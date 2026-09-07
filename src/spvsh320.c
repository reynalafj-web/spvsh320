/*
 * SpVSH-320-FR2 — portable C99 implementation.
 * Matches the unseeded core in the manuscript C++ / ESP32 / SMHasher3 adapters.
 */

#include "spvsh320.h"
#include <string.h>

static const uint64_t SPVSH_P1 = UINT64_C(0xFFFFFFFFFFFFFFC5);
static const uint64_t SPVSH_C1 = 59;
static const uint64_t SPVSH_P2 = UINT64_C(0xFFFFFFFFFFFFFFAD);
static const uint64_t SPVSH_C2 = 83;

static const uint16_t SPVSH_SMALL_PRIMES[64] = {
      3,   5,   7,  11,  13,  17,  19,  23,  29,  31,  37,  41,  43,  47,  53,  59,
     61,  67,  71,  73,  79,  83,  89,  97, 101, 103, 107, 109, 113, 127, 131, 137,
    139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227,
    229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313
};

static int spvsh_ctz64(uint64_t x)
{
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_ctzll(x);
#else
    int n = 0;
    if ((x & 0xFFFFFFFFULL) == 0) { n += 32; x >>= 32; }
    if ((x & 0xFFFFULL) == 0)     { n += 16; x >>= 16; }
    if ((x & 0xFFULL) == 0)       { n += 8;  x >>= 8;  }
    if ((x & 0xFULL) == 0)        { n += 4;  x >>= 4;  }
    if ((x & 0x3ULL) == 0)        { n += 2;  x >>= 2;  }
    if ((x & 0x1ULL) == 0)        { n += 1; }
    return n;
#endif
}

static void mul64x64_128(uint64_t a, uint64_t b, uint64_t *hi, uint64_t *lo)
{
    uint32_t a_lo = (uint32_t)a;
    uint32_t a_hi = (uint32_t)(a >> 32);
    uint32_t b_lo = (uint32_t)b;
    uint32_t b_hi = (uint32_t)(b >> 32);

    uint64_t t = (uint64_t)a_lo * b_lo;
    uint32_t w3 = (uint32_t)t;
    uint32_t k  = (uint32_t)(t >> 32);

    t = (uint64_t)a_hi * b_lo + k;
    uint32_t w2 = (uint32_t)t;
    uint32_t w1 = (uint32_t)(t >> 32);

    t = (uint64_t)a_lo * b_hi + w2;
    k = (uint32_t)(t >> 32);

    *lo = ((uint64_t)(uint32_t)t << 32) | w3;
    *hi = (uint64_t)a_hi * b_hi + w1 + k;
}

static uint64_t mul_mod_p(uint64_t a, uint64_t b, uint64_t P, uint64_t C)
{
    uint64_t hi, lo;
    mul64x64_128(a, b, &hi, &lo);

    uint64_t h_lo = (uint32_t)hi;
    uint64_t h_hi = hi >> 32;
    uint64_t p1 = h_lo * C;
    uint64_t p2 = h_hi * C;
    uint64_t l2 = p1 + (p2 << 32);
    uint64_t h2 = (p2 >> 32) + (l2 < p1 ? 1 : 0);
    uint64_t sum_lo = l2 + lo;
    uint64_t carry  = (sum_lo < lo) ? 1 : 0;
    uint64_t sum_hi = h2 + carry;
    uint64_t final_res = sum_lo + sum_hi * C;
    if (final_res < sum_lo) final_res += C;
    if (final_res >= P) final_res -= P;
    return final_res;
}

static uint64_t load_u64_le(const uint8_t *p)
{
    uint64_t v;
    memcpy(&v, p, 8);
    return v;
}

static void store_u64_le(uint64_t v, uint8_t *p)
{
    memcpy(p, &v, 8);
}

static void permute(spvsh320_ctx *ctx)
{
    uint64_t trigger = ctx->s[0] ^ ctx->s[2] ^ ctx->s[4] ^ UINT64_C(0x9E3779B97F4A7C15);
    uint64_t prod_p1 = 1;
    uint64_t prod_p2 = 1;

    while (trigger != 0) {
        int idx = spvsh_ctz64(trigger);
        uint64_t p = (uint64_t)SPVSH_SMALL_PRIMES[idx];
        prod_p1 = mul_mod_p(prod_p1, p, SPVSH_P1, SPVSH_C1);
        prod_p2 = mul_mod_p(prod_p2, p, SPVSH_P2, SPVSH_C2);
        trigger &= trigger - 1;
    }

    ctx->s[0] = mul_mod_p(mul_mod_p(ctx->s[0], ctx->s[0], SPVSH_P1, SPVSH_C1), prod_p1, SPVSH_P1, SPVSH_C1);
    ctx->s[1] = mul_mod_p(mul_mod_p(ctx->s[1], ctx->s[1], SPVSH_P2, SPVSH_C2), prod_p2, SPVSH_P2, SPVSH_C2);
    ctx->s[2] = mul_mod_p(mul_mod_p(ctx->s[2], ctx->s[2], SPVSH_P1, SPVSH_C1), prod_p1, SPVSH_P1, SPVSH_C1);
    ctx->s[3] = mul_mod_p(mul_mod_p(ctx->s[3], ctx->s[3], SPVSH_P2, SPVSH_C2), prod_p2, SPVSH_P2, SPVSH_C2);
    ctx->s[4] = mul_mod_p(mul_mod_p(ctx->s[4], ctx->s[4], SPVSH_P1, SPVSH_C1), prod_p1, SPVSH_P1, SPVSH_C1);

    {
        uint64_t t0 = ctx->s[0], t1 = ctx->s[1], t2 = ctx->s[2], t3 = ctx->s[3], t4 = ctx->s[4];
        ctx->s[0] = t0 ^ (t1 >> 17) ^ (t1 << 47);
        ctx->s[1] = t1 ^ (t2 >> 11) ^ (t2 << 53);
        ctx->s[2] = t2 ^ (t3 >> 19) ^ (t3 << 45);
        ctx->s[3] = t3 ^ (t4 >> 13) ^ (t4 << 51);
        ctx->s[4] = t4 ^ (t0 >> 23) ^ (t0 << 41);
    }
}

void spvsh320_init(spvsh320_ctx *ctx)
{
    uint64_t seed64 = UINT64_C(0x0123456789ABCDEF);
    int i;

    ctx->s[0] = seed64 + UINT64_C(0xDEADBEEFCAFEBABE);
    ctx->s[1] = seed64 ^ UINT64_C(0x13198A2E03707344);
    ctx->s[2] = seed64 + UINT64_C(0x9E3779B97F4A7C15);
    ctx->s[3] = seed64 ^ UINT64_C(0xBF58476D1CE4E5B9);
    ctx->s[4] = seed64 + UINT64_C(0x94D049BB133111EB);

    for (i = 0; i < SPVSH320_STATE_LANES; ++i) {
        if (ctx->s[i] == 0) ctx->s[i] = 3;
    }
    memset(ctx->buf, 0, sizeof ctx->buf);
    ctx->buf_len = 0;
    permute(ctx);
}

static void absorb_u64(spvsh320_ctx *ctx, uint64_t block)
{
    ctx->s[0] ^= block;
    permute(ctx);
}

void spvsh320_update(spvsh320_ctx *ctx, const uint8_t *data, size_t len)
{
    if (ctx->buf_len > 0) {
        size_t need = SPVSH320_RATE_BYTES - ctx->buf_len;
        size_t take = len < need ? len : need;
        memcpy(ctx->buf + ctx->buf_len, data, take);
        ctx->buf_len += take;
        data += take;
        len  -= take;
        if (ctx->buf_len == SPVSH320_RATE_BYTES) {
            absorb_u64(ctx, load_u64_le(ctx->buf));
            ctx->buf_len = 0;
        }
    }
    while (len >= SPVSH320_RATE_BYTES) {
        absorb_u64(ctx, load_u64_le(data));
        data += SPVSH320_RATE_BYTES;
        len  -= SPVSH320_RATE_BYTES;
    }
    if (len > 0) {
        memcpy(ctx->buf, data, len);
        ctx->buf_len = len;
    }
}

void spvsh320_final(spvsh320_ctx *ctx, uint8_t out[SPVSH320_DIGEST_BYTES])
{
    uint8_t tail[8];
    unsigned r;

    memset(tail, 0, sizeof tail);
    if (ctx->buf_len > 0) memcpy(tail, ctx->buf, ctx->buf_len);
    tail[ctx->buf_len] = 0x80;
    absorb_u64(ctx, load_u64_le(tail));

    for (r = 0; r < SPVSH320_FINALIZATION_ROUNDS; ++r) permute(ctx);

    store_u64_le(ctx->s[0], out + 0);
    permute(ctx);
    store_u64_le(ctx->s[0], out + 8);
    permute(ctx);
    store_u64_le(ctx->s[0], out + 16);
    permute(ctx);
    store_u64_le(ctx->s[0], out + 24);
}

void spvsh320_hash(const uint8_t *data, size_t len, uint8_t out[SPVSH320_DIGEST_BYTES])
{
    spvsh320_ctx ctx;
    spvsh320_init(&ctx);
    spvsh320_update(&ctx, data, len);
    spvsh320_final(&ctx, out);
}

uint64_t spvsh320_first64_be(const uint8_t digest[SPVSH320_DIGEST_BYTES])
{
    uint64_t v = 0;
    int i;
    for (i = 0; i < 8; ++i) v = (v << 8) | digest[i];
    return v;
}
