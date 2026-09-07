# Files that are not in this snapshot yet

This folder currently archives truncated preimage, truncated collision, NIST STS, SMHasher3, ESP32 sketches, and FNIRSI FNB58 traces. If the submitted PDF also cites the items below, add them before you freeze the Zenodo record, or publish them as `paper-v1.1`.

| Cited in the manuscript | Status in this zip | What to add |
|---|---|---|
| Reduced-target preimage search | Present | nothing |
| Reduced-target collision / Pollard's rho | Present | nothing |
| NIST SP 800-22 | Reports + generators present; 1 Gbit streams not stored | optional data-only deposit of `.bin` files |
| SMHasher3 | Adapters + full logs present | optional SMHasher3 commit hash |
| ESP32 throughput and memory | Sketches + Serial logs present | optional `BOARD.md` photo / exact DevKit revision |
| FNIRSI FNB58 energy | `.cfn` + converted CSV present | optional official Toolbox screenshot of Energy_Wh |
| Known-answer tests / standalone library | Present | optional extra vectors beyond tests/test_vectors.txt |
| Figures `alfaj1.png` … | Missing from this software repo | keep figures with the paper, not required here |

Do not invent placeholder CSVs. An incomplete but honest snapshot is better than dummy numbers.
