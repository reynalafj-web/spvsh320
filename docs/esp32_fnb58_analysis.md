# ESP32 latency and FNIRSI FNB58 board-level energy

Board: ESP32 DevKit, 240 MHz, Wi-Fi and Bluetooth off. Sketches compile in Arduino IDE. Heap delta is 0 on every size; the implementations do not allocate.

FNB58 measures USB input of the whole board. It does not isolate the hash core.

---

## Latency (Serial Monitor CSV)

Common sizes only. VSH-256 used 1–10 iterations, so those rows are noisier than the 250–2000 iteration runs.

| Hash | 16 B avg µs | 16 B KB/s | 128 B avg µs | 128 B KB/s | 242 B avg µs | 1024 B avg µs | 1024 B KB/s | 1024 B cpb |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| VSH-256 | 1186 | 13.2 | 8584 | 14.6 | 16 162 | 66 874 | 15.0 | 15 674 |
| SpVSH-320 | 154 | 101 | 392 | 319 | 632 | 2297 | 435 | 538 |
| SHA-3-256 | 258 | 60.5 | 258 | 484 | 515 | 2059 | 486 | 483 |

SHA-3-256 stays near 258 µs for every message that fits in one 136-byte rate block (16, 64, 128 B), then steps when a second absorb starts. SpVSH-320 permutes every 8-byte block, so time grows more smoothly. VSH-256 is a portable double-and-add 256-bit modular multiply; it is the intended slow baseline, not an optimized big-integer port.

At 1024 B, SHA-3-256 is about 1.12× the throughput of SpVSH-320. At 16 B, SpVSH-320 is about 1.67× SHA-3-256.

---

## FNB58 traces

Parsed with the didim99 CFN layout (sample interval 5 s, header field 0.2). Official last-sample Energy_Wh does **not** match ∫V·I dt. Capacity_Ah does match ∫I dt. Report capacity and active-window power; treat the meter Energy_Wh column as a device accumulator that we cannot reconcile.

| Hash | Recorded span (s) | Active span I≥70 mA (s) | I_active (A) | P_active (W) | I_idle (A) | P_idle (W) | P_active − P_idle (W) | Meter last Ah | Meter last Wh | ∫V·I (Wh) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| VSH-256 | 5070 | 3585 | 0.0828 | 0.416 | 0.0538 | 0.270 | 0.146 | 0.1048 | 0.1316 | 0.526 |
| SpVSH-320 | 3635 | 3600 | 0.0818 | 0.411 | 0.0549 | 0.276 | 0.135 | 0.0824 | 0.1034 | 0.413 |
| SHA-3-256 | 3800 | 3570 | 0.0845 | 0.424 | 0.0545 | 0.274 | 0.150 | 0.0874 | 0.1097 | 0.439 |

VSH-256 kept logging ~25 min after the current dropped, so its total watt-hours are not a one-hour figure.

Idle board draw is ~0.27 W. Incremental draw while hashing is only 0.13–0.15 W. Algorithm differences at the USB jack are a few tens of milliwatts on top of the ESP32 baseline. Do not convert these traces into “energy per hash” claims without subtracting idle and stating the 5 s sample interval.

---

## Suggested Results sentences

> On an ESP32 DevKit at 240 MHz we hashed synthetic payloads with Wi-Fi and Bluetooth off. VSH-256, implemented with portable 256-bit double-and-add modular multiplication, ran at about 15 KB/s on 1024-byte messages. SpVSH-320-FR2 reached 435 KB/s and SHA-3-256 reached 486 KB/s on the same size. Short messages reverse the last comparison: 16-byte SpVSH-320 averaged 154 µs against 258 µs for SHA-3-256.
>
> A FNIRSI FNB58 recorded board-level USB voltage and current at 5 s intervals during the one-hour endurance setting. Active-window current was 82–85 mA for all three hashes. Subtracting the 54–55 mA idle floor leaves an incremental 0.13–0.15 W. The meter Energy_Wh register does not match the V·I integral on these files; we therefore report current, power, and capacity rather than a single watt-hour ranking.

---

## Caveats

1. Different iteration counts. Do not treat VSH 1024-byte timing (1 iteration) as equally precise as SHA-3 (1000 iterations).
2. VSH on ESP32 is not the Boost `cpp_int` desktop code. Same modulus and primes; different arithmetic engine.
3. FNB58 is USB-input energy of the DevKit, regulators, and UART, not core joules.
4. Energy_Wh ≠ ∫V·I. Document both. Prefer P_active − P_idle.
5. `.cfn` filenames say “with optimization” for SpVSH and SHA-3. That refers to the `-O3` sketches, not a second algorithm variant.
