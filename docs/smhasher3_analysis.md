# SMHasher3 — three-algorithm comparison

SMHasher3 is a quality-of-mixing suite for general-purpose hashes. A pass is evidence that short keys, sparse keys, and single-bit flips do not produce obvious statistical defects. It is **not** a proof of collision resistance, preimage resistance, or cryptographic security.

Adapters:

- `experiments/smhasher3/src/vsh256.cpp` — same fixed-modulus VSH-256 as the reduced-target programs. Seed is prefixed only when `seed != 0`, so `seed = 0` matches the paper construction.
- `experiments/smhasher3/src/spvsh320.cpp` — SpVSH-320-FR2 core. Seed is mixed into the initial state only when `seed != 0`.
- `experiments/smhasher3/src/sha3.cpp` — SMHasher3 built-in SHA-3-256 (CRYPTOGRAPHIC, official verification `0x79AEFB60`).

SHA-3 ran 237 tests. VSH and SpVSH ran 252 tests. Do not treat the two denominators as the same battery.

---

## Official totals

| Hash | Result | Passed | Seconds |
|---|---|---:|---:|
| VSH-256 | **FAIL** | 125 / 252 | 72 658 |
| SpVSH-320 | **pass** | 252 / 252 | 32 053 |
| SHA-3-256 | **pass** | 237 / 237 | 8 671 |

---

## What failed on VSH-256

The suite marks every Avalanche length and every BIC length as a failure. Representative avalanche maxima:

- 3-byte keys: 100.0%
- 4-byte keys: 99.77%
- 5-byte keys: 76.87%
- 6–192-byte keys: still ~11.5–11.7%

A good 256-bit hash stays near 0.5–1% on this metric. 100% means some input bits leave some output bits unchanged.

The Zeroes keyset is structural, not noise. On 204 800 all-zero keys of increasing length the suite reports **102 399 full 256-bit collisions** (expected ≈ 0). Low-half truncations report 204 799 collisions — almost the entire set. That is consistent with VSH mapping many near-zero messages to a small set of smooth residues.

Other failing families: Sparse, Permutation, Text, TwoBytes, PerlinNoise, Bitflip, and the seed-adapted tests (SeedZeroes, SeedAvalanche, SeedBIC, SeedBitflip). Cyclic fails only the 4-cycle / 3-byte case.

Sanity checks still pass. SMHasher3 can hash and is thread-safe; the construction simply does not mix like a general-purpose hash.

---

## SpVSH-320 and SHA-3-256

Both pass their official totals. Avalanche maxima stay below 1% on every listed key length. BIC coefficients stay near 0.01–0.02. The Zeroes keyset produces zero full-width collisions.

The `-log2(p)` histogram for SpVSH-320 has no mass at 25+. VSH-256 has 3 252 entries in the 25+ bin.

---

## Suggested Results paragraph

> We registered the three constructions in SMHasher3. Seed zero on the VSH-256 and SpVSH-320 adapters matches the unseeded functions used in the reduced-target experiments. VSH-256 failed 127 of 252 tests, including every Avalanche and BIC length and the Zeroes keyset (102 399 observed 256-bit collisions among 204 800 keys; expected approximately 0). SpVSH-320 passed all 252 tests. The built-in SHA-3-256 implementation passed all 237 tests that SMHasher3 scheduled for that hash. SMHasher3 measures mixing quality on short and structured keys. It does not replace a collision or preimage argument.

---

## Caveats for the paper

1. Different test counts (252 vs 237) must stay visible in the table.
2. Seed adapters are harness-only. Do not describe VSH-256 or SpVSH-320 as keyed hashes in the main construction.
3. A SMHasher3 pass does not imply NIST-level cryptographic security, and a SMHasher3 fail does not by itself give a 256-bit collision attack on arbitrary messages.
