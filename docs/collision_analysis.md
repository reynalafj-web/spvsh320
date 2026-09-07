# Truncated collision experiment — three-algorithm comparison

Run `run01`. Method: Pollard's rho + Floyd.  
Transition: \(F_t(x)=\mathrm{Right}_t(H(\mathrm{LE64}(\mathrm{salt}_{t,r})\parallel\mathrm{LE64}(x)))\).  
Expected index: \(\sqrt{\pi/2}\cdot 2^{t/2}\).  
Phase cap: \(2^{32}\). Success: 195/195 per algorithm.

The reported metric is the collision index \(I_{\mathrm{coll}}=\mu+\lambda\), not wall-clock time.

---

## Headline ratios (mean \(I_{\mathrm{coll}} / E[I_{\mathrm{coll}}]\))

| \(t\) | reps | VSH-256 | SpVSH-320 | SHA-3-256 |
|---:|---:|---:|---:|---:|
| 16 | 50 | 1.04 | 0.90 | 0.90 |
| 24 | 50 | 0.99 | 1.06 | 1.09 |
| 32 | 50 | 0.95 | 0.90 | 1.04 |
| 40 | 30 | 1.01 | 1.09 | 1.09 |
| 48 | 10 | 1.00 | 1.04 | 1.01 |
| 56 | 5 | 0.84 | 1.20 | 0.63 |

All three stay near the birthday baseline through 48 bits. The \(t=56\) cells have only five repetitions; SHA-3-256 at 0.63 and SpVSH-320 at 1.20 are compatible with the heavy tail of rho meeting times. Do not read them as a shortcut or as extra strength.

Unlike the preimage experiment, **VSH-256 is not anomalous here**. The rightmost / \(\bmod 2^t\) projection plus the 16-byte salted message removes the high-bit serialization effect that inflated the 8-bit and 16-bit preimage ratios.

---

## Interpretation for the manuscript

These results support a narrow claim only:

> Under a common salted 16-byte rho message and a rightmost truncated projection, none of the three implementations produced an obvious collision shortcut relative to \(\sqrt{\pi/2}\,2^{t/2}\) up to 48 bits. The experiment does not bound full 256-bit collision resistance.

Do not write that SpVSH-320 “matches SHA-3 collision resistance.” The windows are truncated, the message family is one encoding, and \(n=5\) at 56 bits is a pilot.

---

## Suggested Results paragraph

> We searched for truncated collisions with Pollard's rho and Floyd's cycle-finding algorithm. The transition used a domain-separated 16-byte message \(\mathrm{LE64}(\mathrm{salt}_{t,r})\parallel\mathrm{LE64}(x)\) and the rightmost \(t\) digest bits. Every repetition succeeded within the \(2^{32}\) phase cap. Mean collision-index ratios relative to \(\sqrt{\pi/2}\,2^{t/2}\) were 1.04, 0.99, 0.95, 1.01, 1.00, and 0.84 for the evaluated VSH-256 implementation; 0.90, 1.06, 0.90, 1.09, 1.04, and 1.20 for SpVSH-320; and 0.90, 1.09, 1.04, 1.09, 1.01, and 0.63 for SHA-3-256, at \(t=16,24,32,40,48,56\). These values remain consistent with a generic birthday model on the tested windows. They do not establish a collision bound for the full digest.

---

## Repository placement

```
experiments/collision/src/          three C++ programs
experiments/collision/run01/        raw details, summary, reports
results/collision_combined_summary.csv
docs/collision_analysis.md
```
