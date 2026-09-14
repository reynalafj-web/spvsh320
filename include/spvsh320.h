/*
 * SpVSH-320 standalone library
 *
 * 320-bit internal state, 256-bit digest, two output-finalization
 * rounds before the first squeezed lane. Unseeded path used in the
 * PeerJ manuscript experiments.
 *
 * License: MIT
 */
#ifndef SPVSH320_H
#define SPVSH320_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define SPVSH320_DIGEST_BYTES        32
#define SPVSH320_STATE_LANES         5
#define SPVSH320_FINALIZATION_ROUNDS 2
#define SPVSH320_RATE_BYTES          8

typedef struct {
    uint64_t s[SPVSH320_STATE_LANES];
    uint8_t  buf[SPVSH320_RATE_BYTES];
    size_t   buf_len;
} spvsh320_ctx;

void spvsh320_init(spvsh320_ctx *ctx);
void spvsh320_update(spvsh320_ctx *ctx, const uint8_t *data, size_t len);
void spvsh320_final(spvsh320_ctx *ctx, uint8_t out[SPVSH320_DIGEST_BYTES]);

void spvsh320_hash(const uint8_t *data, size_t len, uint8_t out[SPVSH320_DIGEST_BYTES]);

/* First 8 digest bytes as a big-endian integer. Matches the
 * truncated-preimage helper first64_from_output_lane_as_baseline_digest(). */
uint64_t spvsh320_first64_be(const uint8_t digest[SPVSH320_DIGEST_BYTES]);

#ifdef __cplusplus
}
#endif

#endif
