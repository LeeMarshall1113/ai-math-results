#!/bin/bash
# Runs verify.py for one package entirely inside Linux/WSL (proof files stay on the Linux
# filesystem; WSL's /mnt/c is far too slow for tens of GB of LRAT).  usage: run_wsl.sh k9 [verify args]
PKG=$1; shift
HERE="$(cd "$(dirname "$0")/.." && pwd)"
export CADICAL=$HOME/ext/cadical/build/cadical
export GATE_LRAT_CHECK=$HOME/lrchk/bin/gate_lrat_check
export CAKE_LPR=$HOME/lrchk/bin/cake_lpr
export CAKE_HEAP_MB=${CAKE_HEAP_MB:-1000}
export WORK=/tmp/lrw_$PKG
[ "$PKG" = k9 ] && export LEVEL1_EXE=$HOME/lrchk/bin/level1
mkdir -p $WORK
exec python3 $HERE/common/verify.py $HERE/$PKG "$@"
