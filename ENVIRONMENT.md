# Machine and tool versions

Fill the blanks before tagging `paper-v1` if the values are known.
Leave a line as `unknown` rather than guessing.

## Desktop / server (preimage, collision, NIST generators, SMHasher3)

| Item | Value |
|---|---|
| Hostname (optional) | |
| CPU model | |
| Logical cores | |
| RAM | |
| OS | |
| `g++ --version` | |
| Boost Multiprecision | headers only / package version |
| SMHasher3 commit or tag | |
| NIST STS package / version | |
| Worker threads used in run01 | 50 |
| Scheduler (SpVSH preimage) | dynamic work stealing |
| Candidate mode | sequential LE64 counter |
| Finalization rounds | 2 |

## ESP32 (latency, throughput, ROM/RAM, energy)

| Item | Value |
|---|---|
| Board | e.g. ESP32-WROOM-32 DevKitC |
| CPU clock | 240 MHz |
| Arduino-ESP32 core version | |
| `gcc` Xtensa flags | `-O3` / sketch pragma |
| Radios during bench | Wi-Fi and Bluetooth off |
| USB power meter | FNIRSI FNB58 |
| Toolbox CSV sample interval | 5 s (official export) |
| Idle current plateau (approx.) | ~54–55 mA |
| Active current plateau (approx.) | ~82–85 mA |

Board-level USB energy is not core-isolated energy. Do not treat the
FNB58 Energy_Wh column as a formal cryptographic cost metric without
the V·I·t cross-check in `docs/esp32_standardized_addendum.md`.

## How to record versions on a Linux host

```bash
uname -a
lscpu | sed -n '1,20p'
g++ --version | head -1
python3 --version
```
