# Standalone SpVSH-320-FR2 library

Files:

- `include/spvsh320.h`
- `src/spvsh320.c`
- `src/spvsh320.py` (reference, not for timing)
- `tests/test_spvsh320.c`
- `tests/test_vectors.txt`
- `tests/Makefile`

```bash
gcc -O2 -std=c99 -Iinclude -c src/spvsh320.c
cd tests && make test
python3 src/spvsh320.py          # prints SHA-style hex of hash(b"abc")
```

`hash("abc")` =

```
b93de3c241637be1fee467413f12e749796c8d239130cc281619e2d750afc7b8
```

That prefix matches the ESP32 Serial KAT and the SMHasher3 seed=0 path.
