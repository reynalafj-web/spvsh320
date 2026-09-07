#!/bin/sh
set -e
CXX="${CXX:-g++}"
FLAGS="-O3 -std=c++17 -pthread"

echo "Building VSH-256 (needs Boost headers)..."
$CXX $FLAGS -o vsh256_preimage vsh256_preimage.cpp

echo "Building SpVSH-320..."
$CXX $FLAGS -o spvsh320_preimage spvsh320_preimage.cpp

echo "Building SHA-3-256..."
$CXX $FLAGS -o sha3_256_preimage sha3_256_preimage.cpp

echo "Done."
ls -l vsh256_preimage spvsh320_preimage sha3_256_preimage
