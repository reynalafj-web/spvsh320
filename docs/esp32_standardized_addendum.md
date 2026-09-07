# Standardized ESP32 harness + official FNB58 CSV

These files replace the earlier short-iteration latency logs for the paper tables.

Sources:

- `experiments/esp32/src/*_Energy_1H_FNIRSI.ino` — common harness, radios off, 240 MHz
- `experiments/esp32/logs/Result_Throughput.txt` — size sweep, min 1 s or 1000 hashes
- `experiments/esp32/logs/Result_RAM_ROM.txt` — sketch size and heap at boot
- `experiments/fnb58/official_csv/` — FNIRSI Toolbox export, 721 rows, 0–3600 s, 5 s step

`BENCH_AVG_POWER_MW` was left at 0.0 in the flashed sketches. Do not use the energy_per_hash columns from Serial (they are zero). Energy comes from the FNB58 CSV.

---

## Throughput (standardized)

| Hash | 16 B µs | 64 B µs | 128 B µs | 256 B µs | 512 B µs | 1024 B µs | 1024 B KB/s | 1024 B cpb |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| VSH-256 | 1223 | 4367 | 8545 | 16 886 | 33 599 | 66 588 | 15.02 | 15 607 |
| SpVSH-320 | 155 | 257 | 393 | 663 | 1206 | 2293 | 436.2 | 537 |
| SHA-3-256 | 259 | 259 | 259 | 516 | 1030 | 2059 | 485.7 | 483 |

SHA3-256("abc") prefix on this board is `3a985da74fe225b2`, the standard SHA-3-256 digest. SpVSH-320("abc") prefix is `b93de3c241637be1`. VSH-256("abc") printed `0000000000000000` on ESP32. Other VSH sizes in the same run have nonzero prefixes. Treat the ESP32 VSH port as a timing baseline, not as a bit-identical copy of the desktop Boost code, until a known-answer file is checked.

---

## ROM / RAM

Same flash chip (4 MiB) and free heap at start on all three flashes.

| Hash | Sketch bytes | Free heap | Min free heap | Max alloc |
|---|---:|---:|---:|---:|
| VSH-256 | 898 224 | 290 284 | 284 372 | 110 580 |
| SpVSH-320 | 899 280 | 290 284 | 284 372 | 110 580 |
| SHA-3-256 | 899 408 | 290 284 | 284 372 | 110 580 |

The 1 KB sketch-size gap is the algorithm object code plus harness constants. Heap does not move during hashing.

---

## Official FNB58 CSV, first 3600 s

721 samples. Active window = current ≥ 70 mA.

| Hash | Last active t (s) | I_active (mA) | P_active (W) | Meter E at last active (Wh) | Meter C at last active (Ah) |
|---|---:|---:|---:|---:|---:|
| VSH-256 | 3585 | 82.8 | 0.416 | 0.1037 | 0.0826 |
| SpVSH-320 | 3600 | 81.8 | 0.411 | 0.1028 | 0.0818 |
| SHA-3-256 | 3570 | 84.5 | 0.424 | 0.1054 | 0.0839 |

Idle floor after the drop is 55 mA / 0.275 W. Incremental hashing power remains 0.13–0.15 W. Meter Energy_Wh still does not equal ∫V·I; Capacity_Ah does equal ∫I. Quote I_active, P_active, and Capacity. Do not rank the three hashes by the 0.103–0.105 Wh register.

---

## Paper sentence

> A common ESP32 harness at 240 MHz, with radios off, hashed 16-byte to 1024-byte payloads. VSH-256 ran at 15 KB/s on 1024-byte messages. SpVSH-320 reached 436 KB/s and SHA-3-256 reached 486 KB/s. Sketch size differed by about 1 KB. Flash and heap headroom were the same. Board-level USB current during the one-hour FNB58 traces stayed between 82 mA and 85 mA for all three hashes, about 28 mA above the idle floor.
