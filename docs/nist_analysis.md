# NIST SP 800-22 / STS run — three-algorithm comparison

Suite: NIST Statistical Test Suite.  
Dataset size claimed by the reports: 1000 sequences × 1 000 000 bits.  
Minimum pass rate printed by STS: 980/1000 (non-excursion tests).  
Uniformity of p-values fails when the reported P-VALUE is marked `*`.

The 1 Gbit binary streams themselves are **not** stored in this snapshot (too large). What is archived: generators, metadata, and the official `finalAnalysisReport_*.txt` files.

---

## Headline

| Test family | VSH-256 | SpVSH-320 | SHA-3-256 |
|---|---|---|---|
| Frequency | **fail** 451/1000 | pass 997/1000 | pass 988/1000 |
| BlockFrequency | fail uniformity | pass 990/1000 | pass 986/1000 |
| CumulativeSums | **fail** ~450–463 | pass | pass |
| Runs | **fail** 799/1000 | pass 995/1000 | pass 994/1000 |
| LongestRun | **fail** 924/1000 | pass 985/1000 | pass 985/1000 |
| Rank | pass 982/1000 | pass 992/1000 | pass 986/1000 |
| FFT | **fail** 0/1000 | pass 985/1000 | pass 989/1000 |
| NonOverlappingTemplate | **148/148 fail** | 2/148 just below 980 | 0 fail |
| OverlappingTemplate | **fail** 485/1000 | pass 988/1000 | pass 995/1000 |
| Universal | pass 987/1000 | pass 985/1000 | pass 988/1000 |
| ApproximateEntropy | **fail** 1/1000 | pass 992/1000 | pass 990/1000 |
| Serial | **fail** 0/1000 | pass | pass |
| LinearComplexity | pass 988/1000 | pass 985/1000 | pass 989/1000 |
| RandomExcursions(+Variant) | pass (262 seq.) | pass (623 seq.) | pass (624 seq.) |
| Failing instances / 188 | **159** | **2** | **0** |

SpVSH-320 and SHA-3-256 behave like a generic 256-bit digest stream under this suite. VSH-256 does not.

The two SpVSH-320 marks are both NonOverlappingTemplate proportions (979/1000 and 978/1000). That is one count below the printed 980 cutoff, with healthy uniformity p-values. Treat them as threshold noise, not a structural defect.

---

## Why VSH-256 collapses

The attached `vsh_final_1Gb.txt` sample is an ASCII bitstream of 256-bit blocks that are almost all leading zeros. That pattern fails Frequency, FFT, Serial, and ApproximateEntropy immediately.

Two generators exist in this deposit and they are **not equivalent**:

1. `generate_vsh_ascii.cpp` — simplified counter × small-prime product, different 256-bit modulus, one modular reduction, no VSH squaring chain. The sample file matches this style.
2. `generate_nist_vsh256_uniform.cpp` — the same VSH-256 used in the preimage and collision programs (fixed \(n=p q\), block squaring, length bits, 32-byte big-endian serialization).

The STS report labels the generator `<vsh_final.txt>`. Until the full 1 Gbit file and the exact command line are frozen, treat the VSH NIST result as evidence that **this particular VSH bitstream is non-random**, not as a verdict on every VSH parameterization. If the paper claims NIST results for the *same* VSH-256 as the attack experiments, regenerate the stream with `generate_nist_vsh256_uniform.cpp` and rerun STS.

SpVSH-320 and SHA-3-256 reports match the domain-separated generators:

- message \(M_{s,i}=\mathrm{LE64}(\mathrm{salt})\parallel\mathrm{LE64}(s)\parallel\mathrm{LE64}(i)\)
- 256-bit digest concatenation
- SpVSH finalization rounds = 2 (`spvsh320_nist_1000x1M.bin.meta.txt`)

---

## Suggested Results paragraph

> We fed NIST SP 800-22 with 1000 sequences of \(10^6\) bits obtained by concatenating 256-bit digests. The evaluated VSH-256 stream failed 159 of 188 reported test instances, including Frequency (451/1000), FFT (0/1000), Serial (0/1000), ApproximateEntropy (1/1000), and every NonOverlappingTemplate instance. SpVSH-320 failed 2 of 188 instances, both isolated NonOverlappingTemplate proportions of 978/1000 and 979/1000, one count below the suite cutoff of 980/1000. SHA-3-256 failed none of the 188 instances. These outcomes are statistical diagnostics of the generated bitstreams. They are not a proof of collision or preimage resistance.

---

## What is not in this snapshot

- The 1 Gbit `.bin` / `.txt` streams (about 125 MB binary each). Keep them on the experiment machine or on a data-only Zenodo record.
- The truncated 8.4 MB `vsh_final_1Gb.txt` attachment is a sample, not the STS input of record. It is omitted on purpose.
