#!/bin/sh
# Builds both checkers on Linux / WSL / macOS (x64 for cake_lpr).
#  1. gate_lrat_check (our own RUP-only LRAT checker, ~140 lines of C++)
#  2. cake_lpr, the CakeML-verified LRAT/LPR checker (Tan, Heule, Myreen), from the
#     pinned upstream commit; the verified part is cake_lpr.S, whose SHA-256 is checked
#     against the value below (= upstream cake_lpr.sha256 at that commit).
set -e
cd "$(dirname "$0")"
mkdir -p bin
${CXX:-g++} -O2 -std=c++17 -o bin/gate_lrat_check gate_lrat_check.cpp
${CXX:-g++} -O2 -std=c++17 -o bin/level1 level1.cpp
CAKE_COMMIT=a36874a8b750b43fe4b385b8ddbf5b033e46a3fa
CAKE_S_SHA=2f3af32d55083839b3fa0e693afd817679c0b8944bef41def05a8b0ec72b7d4a
if [ ! -d cake_lpr_src ]; then
  git clone -q https://github.com/tanyongkiam/cake_lpr.git cake_lpr_src
fi
( cd cake_lpr_src && git -c advice.detachedHead=false checkout -q $CAKE_COMMIT )
echo "$CAKE_S_SHA  cake_lpr_src/cake_lpr.S" | sha256sum -c -
${CC:-gcc} -O2 -std=c99 cake_lpr_src/basis_ffi.c cake_lpr_src/cake_lpr.S -o bin/cake_lpr 2>/dev/null
bin/cake_lpr --CML_HEAP_SIZE=200 --CML_STACK_SIZE=100 cake_lpr_src/example.cnf cake_lpr_src/example.lpr
echo built: bin/gate_lrat_check bin/cake_lpr
