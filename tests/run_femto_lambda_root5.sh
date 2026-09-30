#!/bin/bash
# Local validation only. Explicit mainconf; no setup.sh/.current_mainconf mutation.
set -euo pipefail
if [[ $# != 2 ]]; then
  echo 'Usage: bash tests/run_femto_lambda_root5.sh MAINCONF ROOT_MACRO_INVOCATION' >&2
  exit 2
fi
test_project=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)
cd "$test_project"
test_tag=$(python3 script/analysis_info_helper.py --library-tag --mainconf "$1")
if [[ "$test_tag" != SL24y ]]; then
  echo 'This regression runner is validated for SL24y/ROOT 5 only.' >&2
  exit 2
fi
test_root=/cvmfs/star.sdcc.bnl.gov/star-spack/spack/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/root-5.34.38-eyawaunrtayjlt3aqon5claakm3c7xl7
exec singularity exec --pwd "$test_project" \
  -B /gpfs:/gpfs -B /cvmfs:/cvmfs -B /star/nfs4/AFS:/star/nfs4/AFS \
  -B /star/data22:/star/data22 -B /star/data24:/star/data24 \
  --env LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}" --env ROOTSYS="$test_root" --env STAR=/star/nfs4/AFS/star/packages/SL24y \
  --env STAR_HOST_SYS=sl73_x8664_gcc485 \
  /cvmfs/singularity.opensciencegrid.org/star-bnl/star-sw:latest /bin/bash -lc '
    export STAR_LIB="$STAR/.$STAR_HOST_SYS/lib"
    export PATH="$ROOTSYS/bin:$STAR/.$STAR_HOST_SYS/bin:$STAR/mgr:$PATH"
    export LD_LIBRARY_PATH="$PWD/lib:$STAR_LIB:$ROOTSYS/lib:${LD_LIBRARY_PATH:-}"
    exec root4star -b -q "$1"
  ' femto-lambda-regression "$2"
