"""SpVSH-320-FR2 reference in Python 3.

Matches src/spvsh320.c and the unseeded C++ / ESP32 / SMHasher3 adapters.
Not a timing-safe or production implementation.
"""
from __future__ import annotations

P1 = (1 << 64) - 59
C1 = 59
P2 = (1 << 64) - 83
C2 = 83
MASK64 = (1 << 64) - 1
GOLDEN = 0x9E3779B97F4A7C15
SMALL_PRIMES = [
      3,   5,   7,  11,  13,  17,  19,  23,  29,  31,  37,  41,  43,  47,  53,  59,
     61,  67,  71,  73,  79,  83,  89,  97, 101, 103, 107, 109, 113, 127, 131, 137,
    139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227,
    229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313,
]
FINALIZATION_ROUNDS = 2


def _mul_mod(a: int, b: int, P: int, C: int) -> int:
    prod = a * b
    lo = prod & MASK64
    hi = prod >> 64
    h_lo = hi & 0xFFFFFFFF
    h_hi = hi >> 32
    p1 = h_lo * C
    p2 = h_hi * C
    l2 = p1 + ((p2 & 0xFFFFFFFF) << 32)
    h2 = (p2 >> 32) + (1 if (l2 & MASK64) < (p1 & MASK64) else 0)
    l2 &= MASK64
    sum_lo = (l2 + lo) & MASK64
    carry = 1 if sum_lo < lo else 0
    sum_hi = h2 + carry
    final_res = (sum_lo + sum_hi * C) & MASK64
    if final_res < sum_lo:
        final_res = (final_res + C) & MASK64
    if final_res >= P:
        final_res -= P
    return final_res


def _permute(s: list[int]) -> None:
    trigger = s[0] ^ s[2] ^ s[4] ^ GOLDEN
    prod_p1 = 1
    prod_p2 = 1
    while trigger:
        idx = (trigger & -trigger).bit_length() - 1
        p = SMALL_PRIMES[idx]
        prod_p1 = _mul_mod(prod_p1, p, P1, C1)
        prod_p2 = _mul_mod(prod_p2, p, P2, C2)
        trigger &= trigger - 1
    s[0] = _mul_mod(_mul_mod(s[0], s[0], P1, C1), prod_p1, P1, C1)
    s[1] = _mul_mod(_mul_mod(s[1], s[1], P2, C2), prod_p2, P2, C2)
    s[2] = _mul_mod(_mul_mod(s[2], s[2], P1, C1), prod_p1, P1, C1)
    s[3] = _mul_mod(_mul_mod(s[3], s[3], P2, C2), prod_p2, P2, C2)
    s[4] = _mul_mod(_mul_mod(s[4], s[4], P1, C1), prod_p1, P1, C1)
    t0, t1, t2, t3, t4 = s
    s[0] = (t0 ^ (t1 >> 17) ^ ((t1 << 47) & MASK64)) & MASK64
    s[1] = (t1 ^ (t2 >> 11) ^ ((t2 << 53) & MASK64)) & MASK64
    s[2] = (t2 ^ (t3 >> 19) ^ ((t3 << 45) & MASK64)) & MASK64
    s[3] = (t3 ^ (t4 >> 13) ^ ((t4 << 51) & MASK64)) & MASK64
    s[4] = (t4 ^ (t0 >> 23) ^ ((t0 << 41) & MASK64)) & MASK64


def _init() -> list[int]:
    seed = 0x0123456789ABCDEF
    s = [
        (seed + 0xDEADBEEFCAFEBABE) & MASK64,
        seed ^ 0x13198A2E03707344,
        (seed + 0x9E3779B97F4A7C15) & MASK64,
        seed ^ 0xBF58476D1CE4E5B9,
        (seed + 0x94D049BB133111EB) & MASK64,
    ]
    s = [3 if x == 0 else x for x in s]
    _permute(s)
    return s


def hash(data: bytes) -> bytes:
    s = _init()
    i = 0
    while i + 8 <= len(data):
        block = int.from_bytes(data[i:i + 8], "little")
        s[0] ^= block
        _permute(s)
        i += 8
    tail = bytearray(8)
    rem = data[i:]
    tail[:len(rem)] = rem
    tail[len(rem)] = 0x80
    s[0] ^= int.from_bytes(tail, "little")
    _permute(s)
    for _ in range(FINALIZATION_ROUNDS):
        _permute(s)
    out = bytearray()
    for k in range(4):
        out += (s[0] & MASK64).to_bytes(8, "little")
        if k < 3:
            _permute(s)
    return bytes(out)


if __name__ == "__main__":
    print(hash(b"abc").hex())
