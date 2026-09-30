#!/bin/bash
# Sourced helper: explicit mainconf; no configuration activation or host starver side effects.
femto_lambda_exec_root() {
  local mode="$1" project="$2" mainconf="$3" root_call="$4"
  shift 4
  local tag star_dir host_sys="" candidate root_sys image path dir seen
  tag="$(python3 "$project/script/analysis_info_helper.py" --library-tag --mainconf "$mainconf")" || return 2
  [[ "$tag" =~ ^[A-Za-z0-9_.-]+$ ]] || { echo "ERROR: invalid libraryTag" >&2; return 2; }
  star_dir="/star/nfs4/AFS/star/packages/$tag"
  for candidate in sl73_x8664_gcc485 sl74_x8664_gcc485; do
    if [[ -d "$star_dir/.$candidate" ]]; then host_sys="$candidate"; break; fi
  done
  [[ -n "$host_sys" ]] || { echo "ERROR: no STAR GCC4.8.5 64-bit build for $tag" >&2; return 2; }
  local roots=()
  shopt -s nullglob
  roots=(/cvmfs/star.sdcc.bnl.gov/star-spack/spack/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/root-5.34.38-*/bin/root-config)
  shopt -u nullglob
  [[ ${#roots[@]} -eq 1 ]] || { echo "ERROR: cannot uniquely resolve batch-matched ROOT5.34.38" >&2; return 2; }
  root_sys="${roots[0]%/bin/root-config}"
  echo "[runtime] STAR=$star_dir STAR_HOST_SYS=$host_sys ROOTSYS=$root_sys"
  if [[ "$mode" == direct ]]; then
    export STAR="$star_dir" STAR_HOST_SYS="$host_sys" ROOTSYS="$root_sys"
    export STAR_LIB="$STAR/.$STAR_HOST_SYS/lib" STAR_BIN="$STAR/.$STAR_HOST_SYS/bin"
    export PATH="$ROOTSYS/bin:$STAR_BIN:$STAR/mgr:$PATH"
    export LD_LIBRARY_PATH="$project/lib:$STAR_LIB:$ROOTSYS/lib:${LD_LIBRARY_PATH:-}"
    exec root4star -b -q "$root_call"
  fi
  command -v singularity >/dev/null || { echo "ERROR: singularity not found" >&2; return 2; }
  image="/cvmfs/singularity.opensciencegrid.org/star-bnl/star-sw:latest"
  local binds=(-B /gpfs:/gpfs -B /cvmfs:/cvmfs -B /star/nfs4/AFS:/star/nfs4/AFS)
  seen="|/gpfs|/cvmfs|/star/nfs4/AFS|"
  for path in "$@"; do
    case "$path" in
      /star/data*/*) dir="/$(printf '%s' "$path" | cut -d/ -f2-3)";;
      /direct/*) dir=/direct;;
      /home/starlib/*) dir=/home/starlib;;
      /tmp/*) dir=/tmp;;
      *) continue;;
    esac
    case "$seen" in *"|$dir|"*) continue;; esac
    [[ -d "$dir" ]] || { echo "ERROR: missing bind directory $dir" >&2; return 2; }
    binds+=(-B "$dir:$dir"); seen="${seen}${dir}|"
  done
  exec singularity exec --pwd "$project" "${binds[@]}" \
    --env "ROOTSYS=$root_sys" --env "STAR=$star_dir" --env "STAR_HOST_SYS=$host_sys" \
    --env "LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-}" \
    "$image" /bin/bash -lc '
      set -euo pipefail
      export STAR_LIB="$STAR/.$STAR_HOST_SYS/lib" STAR_BIN="$STAR/.$STAR_HOST_SYS/bin"
      export PATH="$ROOTSYS/bin:$STAR_BIN:$STAR/mgr:$PATH"
      export LD_LIBRARY_PATH="$1/lib:$STAR_LIB:$ROOTSYS/lib:${LD_LIBRARY_PATH:-}"
      cd "$1"
      exec root4star -b -q "$2"
    ' femto-lambda "$project" "$root_call"
}
