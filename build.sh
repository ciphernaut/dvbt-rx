#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Build dvbt-rx against a checkout of SDRangel (https://github.com/f4exb/sdrangel).
#
#   ./build.sh <sdrangel checkout> [output file]
#
# Needs g++ (C++17) and FFTW's single-precision library (Debian: libfftw3-dev).
# Only SDRangel's modemdvb sources for DVB-T reception are compiled; no Qt.
set -euo pipefail
src="${1:?usage: build.sh <sdrangel checkout> [output file]}"
here="$(cd "$(dirname "$0")" && pwd)"
out="${2:-$here/dvbt-rx}"
m="$src/modemdvb"
[[ -f "$m/dvbtofdm.cpp" ]] || { echo "build.sh: no modemdvb/dvbtofdm.cpp under $src" >&2; exit 1; }
obj="$out.obj"
rm -rf "$obj"
mkdir -p "$obj"
trap 'rm -rf "$obj"' EXIT
flags=(-O2 -std=c++17 -Dsdrangel_STATIC "-I$here/shim" "-I$m" "-I$src/exports")
x86=()
case "$(uname -m)" in
    x86_64|i?86) flags+=(-DMODEMDVB_X86); x86=(dvbtofdm_avx2 resampler_avx2) ;;
esac
for f in dvbt dvbtofdm dvbtinner cpufeatures resampler cwcanceller; do
    g++ "${flags[@]}" -c "$m/$f.cpp" -o "$obj/$f.o"
done
# Compiled for AVX2 and FMA, and chosen at run time when the CPU has them.
for f in "${x86[@]}"; do
    g++ "${flags[@]}" -mavx2 -mfma -c "$m/$f.cpp" -o "$obj/$f.o"
done
g++ "${flags[@]}" -c "$here/main.cpp" -o "$obj/main.o"
g++ "$obj"/*.o -lfftw3f -lpthread -o "$out"
echo "built $out"
