#!/bin/bash
# Finish the proton TOF-threshold scan: merge each point, build its correlation function under its
# own analysis name, and leave the sidecars ready to compare.
#
# The scan varies only pTofPThr, over the same 35 merged blocks, so every point sees the same
# events and the comparison between points is exact. The baseline is the nominal (threshold 0.0)
# over those same 35 blocks, which has to be merged separately from the full nominal set.
#
# Every checkHistAnaFemtoPhi run gets a name suffix: without one they all write the same filenames
# and each run silently replaces the last.
set -u
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
ROOTDIR=/gpfs/mnt/gpfs01/star/pwg/oura/rootfile
IMG=/cvmfs/singularity.opensciencegrid.org/star-bnl/star-sw:latest
MAINCONF=config/mainconf/main_auau3p85fxt_anaFemtoPhiTree_prod_cmp40.yaml
cd "$REPO" || exit 1

hadd_into() {   # hadd_into OUT LIST
    singularity exec -B /gpfs:/gpfs -B /star/nfs4/AFS:/star/nfs4/AFS "$IMG" /bin/bash -lc \
        "export LD_LIBRARY_PATH=\${ROOTSYS}/lib:\$LD_LIBRARY_PATH; cd /tmp; hadd -f $1 \$(cat $2)" \
        > /tmp/scan_hadd.log 2>&1
}

echo "=== $(date +%H:%M:%S) waiting for the scan jobs ==="
while true; do
    n=$(condor_q -totals 2>&1 | grep -o 'oura@bnl.gov: [0-9]*' | awk '{print $2}')
    [ "${n:-1}" = "0" ] && break
    sleep 60
done
echo "=== $(date +%H:%M:%S) jobs done ==="

# Baseline: the same 35 blocks from the nominal run.
NOMDIR=$ROOTDIR/auau3p85fxt_anaFemtoPhiTreeDownstream_prod
find "$NOMDIR" -name 'downstream_block*.root' | sort | head -35 > /tmp/scan_nom35.txt
echo "baseline: $(wc -l < /tmp/scan_nom35.txt) blocks"
hadd_into "$NOMDIR/nominal35_ALL.root" /tmp/scan_nom35.txt
./script/singularity_checkHistAnaFemtoPhi.sh "$NOMDIR/nominal35_ALL.root" "$MAINCONF" _ptof00 \
    > /tmp/scan_chk_00.log 2>&1
echo "  baseline CF: exit $? $(grep -c 'kstarMassFitCF sidecar' /tmp/scan_chk_00.log) sidecar"

for thr in 1.0 1.5 2.0 3.0; do
    tag=$(echo "$thr" | tr -d '.')
    D=$ROOTDIR/auau3p85fxt_anaFemtoPhiTreeDownstream_ptof$tag
    n=$(find "$D" -name '*.root' -size +1M 2>/dev/null | wc -l)
    echo "=== pTofPThr=$thr : $n outputs ==="
    if [ "$n" -eq 0 ]; then echo "  no outputs, skipping"; continue; fi
    find "$D" -name 'downstream_ptof*_block*.root' | sort > /tmp/scan_list_$tag.txt
    hadd_into "$D/ptof${tag}_ALL.root" /tmp/scan_list_$tag.txt
    ./script/singularity_checkHistAnaFemtoPhi.sh "$D/ptof${tag}_ALL.root" "$MAINCONF" _ptof$tag \
        > /tmp/scan_chk_$tag.log 2>&1
    echo "  CF: exit $?  $(grep 'kstarMassFitCF sidecar' /tmp/scan_chk_$tag.log | head -1)"
done
echo "=== $(date +%H:%M:%S) scan summary ready ==="
