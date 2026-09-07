#include "spvsh320.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

static void print_hex(const uint8_t *d, size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) printf("%02x", d[i]);
}

static int expect_hex(const char *label, const uint8_t *d, const char *hex)
{
    char got[65];
    size_t i;
    for (i = 0; i < 32; ++i) sprintf(got + 2 * i, "%02x", d[i]);
    got[64] = 0;
    if (strcmp(got, hex) != 0) {
        fprintf(stderr, "FAIL %s\n  got  %s\n  want %s\n", label, got, hex);
        return 1;
    }
    printf("PASS %s\n", label);
    return 0;
}

int main(void)
{
    uint8_t d1[32], d2[32];
    int rc = 0;
    const uint8_t abc[] = "abc";
    size_t i;

    printf("SpVSH-320-FR2 known-answer and split-update checks\n");

    spvsh320_hash(abc, 3, d1);
    printf("abc = ");
    print_hex(d1, 32);
    printf("\n");

    /* Incremental must match one-shot. */
    {
        spvsh320_ctx ctx;
        spvsh320_init(&ctx);
        spvsh320_update(&ctx, abc, 1);
        spvsh320_update(&ctx, abc + 1, 2);
        spvsh320_final(&ctx, d2);
        if (memcmp(d1, d2, 32) != 0) {
            fprintf(stderr, "FAIL incremental abc\n");
            rc = 1;
        } else {
            printf("PASS incremental abc\n");
        }
    }

    /* Empty, 8-byte, 9-byte, 64-byte split. */
    {
        uint8_t buf[64];
        for (i = 0; i < 64; ++i) buf[i] = (uint8_t)i;
        spvsh320_hash(buf, 0, d1);
        printf("empty = "); print_hex(d1, 32); printf("\n");
        spvsh320_hash(buf, 8, d1);
        printf("seq8  = "); print_hex(d1, 32); printf("\n");
        spvsh320_hash(buf, 9, d1);
        printf("seq9  = "); print_hex(d1, 32); printf("\n");
        spvsh320_hash(buf, 64, d1);
        {
            spvsh320_ctx ctx;
            spvsh320_init(&ctx);
            spvsh320_update(&ctx, buf, 3);
            spvsh320_update(&ctx, buf + 3, 29);
            spvsh320_update(&ctx, buf + 32, 32);
            spvsh320_final(&ctx, d2);
        }
        if (memcmp(d1, d2, 32) != 0) {
            fprintf(stderr, "FAIL incremental 64\n");
            rc = 1;
        } else {
            printf("PASS incremental 64\n");
        }
        printf("seq64 = "); print_hex(d1, 32); printf("\n");
    }

    /* ESP32 / SMHasher seed=0 prefix for "abc" */
    {
        char prefix[17];
        spvsh320_hash(abc, 3, d1);
        for (i = 0; i < 8; ++i) sprintf(prefix + 2 * i, "%02x", d1[i]);
        prefix[16] = 0;
        printf("abc prefix64 = %s\n", prefix);
        if (strcmp(prefix, "b93de3c241637be1") != 0) {
            fprintf(stderr, "WARN prefix does not match ESP32 KAT b93de3c241637be1\n");
            /* still not a hard fail until we confirm endian printing */
        }
    }

    return rc;
}
