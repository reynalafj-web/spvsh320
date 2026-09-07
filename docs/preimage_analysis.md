# Truncated preimage experiment — three-algorithm comparison

Run `run01`, 50 workers, cap \(2^{40}\) trials per repetition.
Target: leftmost \(t\) bits of \(H(\texttt{TRUNCATED\_PREIMAGE\_TARGET\_001|rep=}r)\).
Candidates: 8-byte little-endian counter \(M_i=\mathrm{LE64}(i)\).
Repetitions: 50 for \(t\in\{8,16,24,32\}\), 30 for \(t=36\).
Success: 230/230 for every algorithm (no timeout).

Metric reported in the paper should remain **work-equivalent trial count**, not wall-clock time. Schedulers differ (VSH-256: fixed, one lane per repetition; SpVSH-320 and SHA-3-256: dynamic work stealing), but each reported `trials` value is the number of hashes actually evaluated until the first match.

---

## 1. Headline table (copy into the manuscript)

Mean ratio = mean(observed trials) / \(2^t\). A ratio near 1 is consistent with a generic uniform model. A ratio \(\gg 1\) means the target prefix is harder to hit from the 8-byte counter family than a random function would predict. A ratio \(\ll 1\) would suggest a shortcut.

| Algorithm | \(t\) | Success | Mean trials | Expected \(2^t\) | Mean ratio | Median ratio | Mean \(\log_2\) trials |
|---|---:|---:|---:|---:|---:|---:|---:|
| VSH-256 | 8 | 50/50 | \(1.80\times 10^{5}\) | 256 | **704.72** | 747.88 | 17.40 |
| VSH-256 | 16 | 50/50 | \(1.18\times 10^{6}\) | 65 536 | **18.00** | 15.24 | 20.04 |
| VSH-256 | 24 | 50/50 | \(3.49\times 10^{7}\) | \(1.68\times 10^{7}\) | 2.08 | 1.74 | 24.72 |
| VSH-256 | 32 | 50/50 | \(4.91\times 10^{9}\) | \(4.29\times 10^{9}\) | 1.14 | 0.88 | 31.33 |
| VSH-256 | 36 | 30/30 | \(8.15\times 10^{10}\) | \(6.87\times 10^{10}\) | 1.19 | 0.91 | 35.85 |
| SpVSH-320 | 8 | 50/50 | 268 | 256 | 1.05 | 0.78 | 7.13 |
| SpVSH-320 | 16 | 50/50 | \(5.40\times 10^{4}\) | 65 536 | 0.82 | 0.66 | 14.85 |
| SpVSH-320 | 24 | 50/50 | \(1.56\times 10^{7}\) | \(1.68\times 10^{7}\) | 0.93 | 0.53 | 22.92 |
| SpVSH-320 | 32 | 50/50 | \(4.33\times 10^{9}\) | \(4.29\times 10^{9}\) | 1.01 | 0.86 | 31.24 |
| SpVSH-320 | 36 | 30/30 | \(7.56\times 10^{10}\) | \(6.87\times 10^{10}\) | 1.10 | 0.60 | 35.27 |
| SHA-3-256 | 8 | 50/50 | 304 | 256 | 1.19 | 1.13 | 7.61 |
| SHA-3-256 | 16 | 50/50 | \(7.67\times 10^{4}\) | 65 536 | 1.17 | 0.75 | 15.16 |
| SHA-3-256 | 24 | 50/50 | \(1.90\times 10^{7}\) | \(1.68\times 10^{7}\) | 1.13 | 0.75 | 23.40 |
| SHA-3-256 | 32 | 50/50 | \(4.80\times 10^{9}\) | \(4.29\times 10^{9}\) | 1.12 | 0.77 | 31.14 |
| SHA-3-256 | 36 | 30/30 | \(9.02\times 10^{10}\) | \(6.87\times 10^{10}\) | 1.31 | 0.82 | 35.66 |

---

## 2. What the numbers support

### SpVSH-320
No evidence of a truncated-preimage shortcut in this setting. Mean ratios sit between 0.82 and 1.10. At \(t=32\) and \(t=36\) the cost tracks \(2^t\). The slightly sub-unity means at \(t=16\) and \(t=24\) are consistent with a heavy-tailed geometric sample (\(n=50\)), not with an attack. Do **not** write that SpVSH-320 is “easier than generic.”

### SHA-3-256
Same qualitative picture. Mean ratios 1.12–1.31. The \(t=36\) mean of 1.31 with \(n=30\) is still compatible with geometric variance (the maximum trial count in that cell is \(4.43\times 10^{11}\), about \(6.4\cdot 2^{36}\)). Use SHA-3-256 as the healthy reference, not as a proof that SpVSH matches SHA-3 security.

### VSH-256
This is the result that needs a paragraph, not a footnote.

- At \(t=8\), every one of 50 repetitions needed on the order of \(2^{17}\) trials (mean ratio **705**, min already 383).
- At \(t=16\) the inflation is still large (mean ratio **18**).
- By \(t=24\) it drops to 2.08, and by \(t=32\)–\(36\) it is close to the other two algorithms (1.14–1.19).

So the evaluated VSH-256 implementation does **not** behave like a uniform random function on its **leftmost short prefixes** when the candidate family is 8-byte counters and the target family is the longer ASCII strings used here. The same implementation becomes compatible with the \(2^t\) model once \(t\) is large enough that the structured high bits are no longer the whole truncated window.

This is **not** a proof that VSH-256 is easy to invert. It is the opposite at small \(t\): the target prefix is *harder* to realize from short counters than a random 8-bit string would be. Likely causes, which should be stated as hypotheses:

1. High-bit structure of \(x^2\cdot\prod p_i \bmod n\) for the 256-bit RSA-style modulus used in the implementation (\(k=43\) primes).
2. Length / padding mismatch: targets are ~35-byte ASCII strings; candidates are always 8 bytes. VSH folds bits into prime products blockwise, so different message lengths occupy different regions of the high-bit space.
3. The experiment extracts `first64 = (x >> 192)` and then the leftmost \(t\) bits of that word. Any bias in the top bits of the modular representative is inherited by small \(t\).

Do not hide the \(t=8\) column. Reviewers will ask. Interpret it as a **diagnostic of output regularity**, which is consistent with the weak NIST / SMHasher3 profile already reported for this VSH implementation.

---

## 3. Text you can paste into Results

**Compact paragraph (recommended):**

> All three implementations found a truncated preimage in every repetition at \(t\in\{8,16,24,32,36\}\), with a per-repetition cap of \(2^{40}\) candidate evaluations. For SpVSH-320 the mean work ratio relative to \(2^t\) was 1.05, 0.82, 0.93, 1.01, and 1.10. For SHA-3-256 the corresponding ratios were 1.19, 1.17, 1.13, 1.12, and 1.31. These values remain consistent with a generic model on the tested truncated windows; they do not indicate a shortcut, and they do not establish full-digest preimage resistance. The evaluated VSH-256 implementation behaved differently at short prefixes: the mean ratios were 705 at 8 bits and 18.0 at 16 bits, then fell to 2.08, 1.14, and 1.19 at 24, 32, and 36 bits. The short-prefix inflation is consistent with structured high bits of the modular VSH output under an 8-byte counter candidate family, rather than with a uniform random 8- or 16-bit projection. Because the comparison uses a single candidate encoding and a single target-message family, the ratios should be read as a same-scenario diagnostic and not as a concrete preimage bound.

---

## 4. What to put in the repository

Keep these files together under `experiments/preimage/run01/`:

```
vsh256_preimage.cpp
vsh256_preimage_run01_details.csv
vsh256_preimage_run01_summary.csv
vsh256_preimage_run01_report.txt
spvsh320_preimage.cpp
spvsh320_preimage_run01_details.csv
spvsh320_preimage_run01_summary.csv
spvsh320_preimage_run01_report.txt
sha3_256_preimage.cpp
sha3_256_preimage_run01_details.csv
sha3_256_preimage_run01_summary.csv
sha3_256_preimage_run01_report.txt
preimage_combined_summary.csv
```

Record in the README:

- VSH-256: `--workers 50`, fixed scheduler, `lanes_per_rep=1`, `k=43`, modulus `n=p*q` with the two 128-bit primes in the source.
- SpVSH-320: `--workers 50 --scheduler dynamic --candidate-mode counter --finalization-rounds 2`.
- SHA-3-256: `--workers 50 --scheduler dynamic`, unseeded FIPS 202 path (`seed=0`).
- Host: same class of machine for all three (state the exact CPU if you still have it; the reports say worker budget 50).

---

## 5. Do not claim

- “SpVSH-320 is as preimage-resistant as SHA-3.” The test is truncated, same-length candidates only, and stops at 36 bits.
- “VSH-256 is broken.” Small-\(t\) inflation is a regularity / domain-separation effect, not a cheap inverter of the 256-bit digest.
- “Ratio 0.82 at 16 bits means a 18% shortcut.” Geometric sampling with \(n=50\) easily produces that.
- Comparability of `mean_seconds` across algorithms. VSH uses `boost::cpp_int`; SpVSH and SHA-3 use native 64-bit paths. Seconds measure implementation cost, not security.

---

## 6. Optional extra check (not required before submission)

If a reviewer asks why VSH \(t=8\) is so large, a cheap extra experiment is enough:

- Hash \(2^{20}\) random 8-byte counters and histogram the leftmost 8 bits.
- Hash the 50 ASCII targets and see where those 8-bit prefixes sit in that histogram.

If some target prefixes have frequency \(\approx 2^{-8}/700\), the ratio is explained. You do not need to rerun the \(2^{40}\) search.
