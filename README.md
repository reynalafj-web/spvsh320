# SpVSH-320 reproducibility package

Reference experiment code and raw outputs for:

> Alfajri, R., & Afianti, F. *SpVSH-320: A Sponge-Like Hash Construction with Smooth-Number-Inspired State Transformation.* Manuscript submitted to *PeerJ Computer Science*.

This snapshot is an **experimental research package**. It is not a production cryptographic library.

This release contains the **truncated preimage experiment** that is described in the paper appendix:

- target: leftmost \(t\) bits of \(H(\texttt{TRUNCATED\_PREIMAGE\_TARGET\_001|rep=}r)\)
- candidates: 8-byte little-endian counters \(M_i=\mathrm{LE64}(i)\)
- lengths: \(t \in \{8,16,24,32,36\}\)
- repetitions: 50 for 8–32 bits, 30 for 36 bits
- cap: \(2^{40}\) evaluations per repetition

This snapshot also contains truncated collision results, NIST SP 800-22 reports, SMHasher3 adapters and logs, ESP32 Arduino sketches with Serial timing logs, and FNIRSI FNB58 `.cfn` traces (CSV in `experiments/fnb58/csv/`).

## Build

Requirements:

- C++17 compiler (`g++` or `clang++`)
- Boost.Multiprecision headers for the VSH-256 program only

```bash
sudo apt-get update
sudo apt-get install -y g++ libboost-dev

cd experiments/preimage/src

g++ -O3 -std=c++17 -pthread -o vsh256_preimage vsh256_preimage.cpp
g++ -O3 -std=c++17 -pthread -o spvsh320_preimage spvsh320_preimage.cpp
g++ -O3 -std=c++17 -pthread -o sha3_256_preimage sha3_256_preimage.cpp
```

The paper runs used 50 worker threads on a dedicated server. A laptop can rebuild the binaries and replay the short \(t=8\) cells. The \(t=36\) cells need a long wall-clock budget.

## Reproduce run01

The commands below recreate the configuration recorded in the `*_report.txt` files. They will overwrite CSV names unless you change `--out-prefix`.

```bash
# VSH-256: fixed scheduler, one search lane per repetition
./vsh256_preimage \
  --workers 50 \
  --lanes-per-rep 1 \
  --bits 8,16,24,32,36 \
  --max-trials-power 40 \
  --out-prefix vsh256_preimage_run01

# SpVSH-320: dynamic work stealing, sequential counters, two finalization rounds
./spvsh320_preimage \
  --workers 50 \
  --scheduler dynamic \
  --candidate-mode counter \
  --finalization-rounds 2 \
  --bits 8,16,24,32,36 \
  --max-trials-power 40 \
  --out-prefix spvsh320_preimage_run01

# SHA-3-256: dynamic work stealing, unseeded FIPS 202 path
./sha3_256_preimage \
  --workers 50 \
  --scheduler dynamic \
  --bits 8,16,24,32,36 \
  --max-trials-power 40 \
  --out-prefix sha3_256_preimage_run01
```

The security-oriented metric is the **number of candidate hashes evaluated**, not wall-clock time. VSH-256 uses `boost::cpp_int` and is much slower per trial than the native 64-bit SpVSH-320 and SHA-3-256 paths.

## How to cite

Until the article exists:

```bibtex
@software{spvsh320_repro,
  author  = {Alfajri, Reynaldi and Afianti, Farah},
  title   = {SpVSH-320 reproducibility package},
  year    = {2026},
  version = {paper-v1},
  doi     = {10.5281/zenodo.XXXXXXX},
  url     = {https://doi.org/10.5281/zenodo.XXXXXXX}
}
```

After publication, cite the PeerJ Computer Science article and keep this DOI for the software snapshot.

## Disclaimer

These programs are research instruments. Do not deploy SpVSH-320 as a security dependency. The truncated-preimage ratios are same-scenario diagnostics. They are not a reduction-based bound on the 256-bit digest.
