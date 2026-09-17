#!/bin/bash
# Drive the post-production merge to completion, one wave at a time.
#
# The merged blocks and the subjob outputs they were built from do not both fit under the quota,
# so the work has to be done in waves: merge a range of blocks, verify each against the files it
# was built from, release the inputs of the blocks that passed, and only then start the next wave.
# Each step is a farm submission; this script is only the sequencing between them.
#
# It stops rather than continuing whenever the picture is not clean: an unfinished .part file after
# a merge step, or any block that failed verification. Nothing is released in either case, because
# releasing is the one step in this workflow that cannot be undone without re-running the farm.
#
# Usage: script/run_tree_merge_waves.sh WORKDIR OUTDIR STEM GOODLIST FIRST:LAST[,FIRST:LAST...]
set -u
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
WORKDIR="${1:?workdir}"; OUTDIR="${2:?outdir}"; STEM="${3:?stem}"; GOOD="${4:?good list}"
WAVES="${5:?wave ranges, e.g. 70:105,105:141}"
SUBMIT="python3 $REPO/script/submit_tree_merge.py"

wait_for_queue() {
    # condor_q reports this user's total; the merge and verify steps are the only jobs this
    # workflow submits, so an empty queue means the step has finished one way or another.
    while true; do
        n=$(condor_q -totals 2>&1 | grep -o 'oura@bnl.gov: [0-9]*' | awk '{print $2}')
        [ "${n:-1}" = "0" ] && break
        sleep 60
    done
}

cd "$REPO" || exit 1
for wave in $(echo "$WAVES" | tr ',' ' '); do
    echo "=== $(date +%H:%M:%S) wave $wave: merge ==="
    $SUBMIT merge --good "$GOOD" --outdir "$OUTDIR" --workdir "$WORKDIR" --stem "$STEM" \
        --block 100 --block-range "$wave" || exit 1
    wait_for_queue
    parts=$(find "$OUTDIR" -name '*.part' | wc -l)
    if [ "$parts" != "0" ]; then
        echo "STOP: $parts unfinished .part files after wave $wave"; exit 1
    fi

    echo "=== $(date +%H:%M:%S) wave $wave: verify ==="
    $SUBMIT verify --outdir "$OUTDIR" --workdir "$WORKDIR" --shards 16 || exit 1
    wait_for_queue

    fails=$(cat "$WORKDIR"/verified_*.txt 2>/dev/null | grep -c '^FAIL')
    if [ "$fails" != "0" ]; then
        echo "STOP: $fails blocks failed verification; nothing released"; exit 1
    fi

    echo "=== $(date +%H:%M:%S) wave $wave: release the inputs of verified blocks ==="
    $SUBMIT reclaim --outdir "$OUTDIR" --workdir "$WORKDIR" --apply || exit 1
    df -h /gpfs/mnt/gpfs01 | tail -1
done
echo "=== $(date +%H:%M:%S) all waves done: $(find "$OUTDIR" -name '*.root' | wc -l) blocks ==="
