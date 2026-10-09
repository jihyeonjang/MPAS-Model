#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "$0")/../.." && pwd)"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
cd "$tmp"
"${FC:-gfortran}" -std=f2008 -fcheck=all -Wall -Wextra \
 "$root/src/core_atmosphere/physics/mpas_rrtmgp_thermodynamics.F90" \
 "$root/testing_and_setup/rrtmgp/test_mpas_rrtmgp_thermodynamics.F90" -o test_rrtmgp_conversions
./test_rrtmgp_conversions
