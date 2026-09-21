#!/usr/bin/env python3
"""Submit the reduced-tree downstream reader over every verified block as Condor jobs.

One block is one job. The reader is run through
``script/singularity_run_anaFemtoPhiTreeDownstreamV3.sh`` inside a self-contained runtime bundle
so that editing the working tree while jobs are queued cannot change what they execute.

Two modes:

* ``--mode phi``      the phi-p / phi-d channel set from the mainconf (M(KK) plus the extended
                      pair QA).
* ``--mode data006``  the p-p / p-d / d-d track-track channel set from a DATA-006 config.

The two runs are deliberately given the same pair-mT and pair-rapidity axes, so every QA
distribution can be compared bin by bin between the resonance and the track-track channels.

Typical use::

    # 1. one block, to prove the bundle and the farm environment
    script/submit_downstream_blocks.py --mode data006 --tag data006_full_20260919 --limit 1
    # 2. everything else, once that job has exited 0
    script/submit_downstream_blocks.py --mode data006 --tag data006_full_20260919

The second call reuses the bundle and skips blocks whose output already exists, so it never
resubmits the smoke block.

The skip test is only "an output file is there and is not empty". ROOT creates the output file
when the job starts, so a running, a killed and a finished job all look the same to it. That is
what keeps a smoke job from being resubmitted, but it means the outputs have to be verified
against the job logs -- every log ending in `exit 0`, and the processed-event counts adding up to
the block totals -- before a run is called complete. Use `--force` to redo a slice.
"""

import argparse
import glob
import os
import shlex
import shutil
import subprocess
import sys

PROJECT_ROOT = os.path.realpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
BUNDLE_DIRS = ("config", "include", "lib", "script", "StMaker", "tools")
DEFAULT_VERIFIED_GLOB = os.path.join(PROJECT_ROOT, "job/merge/prod20260916/verified_*.txt")
DEFAULT_MAINCONF = "config/mainconf/main_auau3p85fxt_anaFemtoPhiTree_prod.yaml"
DEFAULT_DATA006_CONFIG = "config/maker/data006_auau3p85fxt_pp_pd_dd.yaml"


def parse_args(argv):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--mode", choices=("phi", "data006"), required=True)
    p.add_argument("--tag", required=True,
                   help="work-directory and output-file tag, e.g. data006_full_20260919")
    p.add_argument("--mainconf", default=DEFAULT_MAINCONF)
    p.add_argument("--data006-config", default=DEFAULT_DATA006_CONFIG,
                   help="only used by --mode data006")
    p.add_argument("--outdir", default=None,
                   help="ROOT output directory (default: "
                        "/gpfs/mnt/gpfs01/star/pwg/oura/rootfile/<tag>)")
    p.add_argument("--workdir", default=None,
                   help="job/log directory (default: job/merge/<tag>)")
    p.add_argument("--verified-glob", default=DEFAULT_VERIFIED_GLOB,
                   help="glob of PASS/FAIL lists that name the input blocks")
    p.add_argument("--block-list", default=None,
                   help="explicit file with one input ROOT path per line, instead of the lists")
    p.add_argument("--max-events", type=int, default=-1)
    p.add_argument("--events-per-job", type=int, default=0,
                   help="split each block into slices of N events, one job per slice (0 = one "
                        "job per whole block). A p-p / p-d / d-d pass over a whole 11 M-event "
                        "block is about a day of wall clock, past what a batch slot will hold.")
    p.add_argument("--pair-mt-bins", type=int, default=150)
    p.add_argument("--pair-mt-min", type=float, default=0.8)
    p.add_argument("--pair-mt-max", type=float, default=3.8)
    p.add_argument("--buffer-size", type=int, default=-1)
    p.add_argument("--mixing-mode", default="bufferAll")
    p.add_argument("--variation", default="")
    p.add_argument("--three-event", action="store_true",
                   help="fill d(E0)-K+(E1)-K-(E2) background; phi mode only")
    p.add_argument("--three-event-cap", type=int,
                   help="maximum raw three-event triplets sampled per current event")
    p.add_argument("--three-event-seed", type=int,
                   help="positive deterministic three-event sampling seed")
    p.add_argument("--request-memory", default="2 GB")
    p.add_argument("--limit", type=int, default=0, help="submit at most N jobs (0 = no limit)")
    p.add_argument("--force", action="store_true",
                   help="also submit blocks whose output already exists")
    p.add_argument("--reuse-runtime", action="store_true",
                   help="keep the existing runtime bundle instead of rebuilding it")
    p.add_argument("--dry-run", action="store_true", help="write the files but do not submit")
    return p.parse_args(argv)


def read_blocks(args):
    if args.block_list:
        with open(args.block_list) as fh:
            paths = [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]
    else:
        paths = []
        for listing in sorted(glob.glob(args.verified_glob)):
            with open(listing) as fh:
                for line in fh:
                    parts = line.split()
                    if len(parts) == 2 and parts[0] == "PASS":
                        paths.append(parts[1])
    missing = [p for p in paths if not os.path.isfile(p)]
    if missing:
        sys.exit("ERROR: %d listed input block(s) do not exist, first: %s"
                 % (len(missing), missing[0]))
    if not paths:
        sys.exit("ERROR: no input blocks found")
    return sorted(set(paths))


SINGULARITY_IMAGE = "/cvmfs/singularity.opensciencegrid.org/star-bnl/star-sw:latest"


def block_event_counts(blocks, cache):
    """{path: event count}, read from `cache` or measured once with the container's own ROOT."""
    counts = {}
    if os.path.isfile(cache):
        with open(cache) as fh:
            for line in fh:
                parts = line.split(None, 1)
                if len(parts) == 2:
                    counts[parts[1].strip()] = int(parts[0])
    todo = [b for b in blocks if b not in counts]
    if todo:
        print("[submit] counting events in %d block(s) ..." % len(todo))
        listing = cache + ".list"
        with open(listing, "w") as fh:
            fh.write("\n".join(todo) + "\n")
        macro = os.path.join(PROJECT_ROOT, "tools/print_tree_event_counts.C")
        cmd = ["singularity", "exec", "-B", "/gpfs:/gpfs", "-B", "/star/u:/star/u",
               SINGULARITY_IMAGE, "/bin/bash", "-lc",
               "root -b -q %s" % shlex.quote("%s(\"%s\")" % (macro, listing))]
        out = subprocess.check_output(cmd, stderr=subprocess.DEVNULL).decode()
        for line in out.splitlines():
            parts = line.split()
            if len(parts) == 3 and parts[0] == "COUNT":
                counts[parts[2]] = int(parts[1])
        os.remove(listing)
        with open(cache, "w") as fh:
            for path in sorted(counts):
                fh.write("%d %s\n" % (counts[path], path))
    bad = [b for b in blocks if counts.get(b, -1) <= 0]
    if bad:
        sys.exit("ERROR: could not read an event count for %d block(s), first: %s"
                 % (len(bad), bad[0]))
    return counts


def block_key(path):
    """block0007 out of ..._block0007.root; falls back to the bare stem."""
    stem = os.path.basename(path)[:-len(".root")] if path.endswith(".root") else \
        os.path.basename(path)
    idx = stem.rfind("block")
    return stem[idx:] if idx >= 0 else stem


def build_runtime(bundle, reuse):
    if reuse and os.path.isdir(bundle):
        print("[submit] reusing runtime bundle %s" % bundle)
        return
    if os.path.isdir(bundle):
        shutil.rmtree(bundle)
    os.makedirs(bundle)
    for name in BUNDLE_DIRS:
        src = os.path.join(PROJECT_ROOT, name)
        if not os.path.isdir(src):
            sys.exit("ERROR: missing bundle directory %s" % src)
        shutil.copytree(src, os.path.join(bundle, name), symlinks=True)
    print("[submit] runtime bundle written to %s" % bundle)


def job_script(args, bundle, tree, out_root, first_event, max_events):
    wrapper = os.path.join(bundle, "script/singularity_run_anaFemtoPhiTreeDownstreamV3.sh")
    data006 = args.data006_config if args.mode == "data006" else ""
    # The wrapper's positional interface; every argument is given explicitly so a job never
    # falls back on a default that the next edit of the wrapper could change underneath it.
    three_enabled = "1" if args.three_event else "0"
    three_cap = str(args.three_event_cap) if args.three_event else "0"
    three_seed = str(args.three_event_seed) if args.three_event else "0"
    call = [wrapper, args.mainconf, tree, out_root,
            str(args.buffer_size), args.mixing_mode, args.variation,
            "-1", "-1", "0", "0", str(max_events),
            str(args.pair_mt_bins), str(args.pair_mt_min), str(args.pair_mt_max),
            "0", data006, str(first_event), three_enabled, three_cap, three_seed]
    quoted = " ".join(shlex.quote(a) for a in call)
    return """#!/bin/bash
set -u
echo "host $(hostname)  start $(date)"
cd %s || exit 1
%s
rc=$?
echo "exit $rc  end $(date)"
if [ $rc -eq 0 ]; then ls -l %s; fi
exit $rc
""" % (bundle, quoted, out_root)


def condor_block(args, workdir, name):
    sh = os.path.join(workdir, "job", "%s.sh" % name)
    return """Universe         = vanilla
Notification     = never
Executable       = /bin/sh
Arguments        = "-c 'exec %s'"
Output           = %s/log/%s.out
Error            = %s/log/%s.err
Log              = %s/log/condor.log
Initialdir       = %s
kill_sig         = SIGINT
request_memory   = %s
PeriodicRemove   = (JobStatus == 2 && (CurrentTime - JobCurrentStartDate > (43200))) || (((CurrentTime - EnteredCurrentStatus) > (2*24*3600)) && JobStatus == 5)
Priority         = +10
Queue

""" % (sh, workdir, name, workdir, name, workdir, workdir, args.request_memory)


def main(argv):
    args = parse_args(argv)
    workdir = os.path.realpath(args.workdir or os.path.join(PROJECT_ROOT, "job/merge", args.tag))
    outdir = os.path.realpath(
        args.outdir or "/gpfs/mnt/gpfs01/star/pwg/oura/rootfile/%s" % args.tag)
    bundle = os.path.join(workdir, "runtime")

    for path in (os.path.join(workdir, "job"), os.path.join(workdir, "log"), outdir):
        if not os.path.isdir(path):
            os.makedirs(path)

    if args.mode == "data006":
        cfg = os.path.join(PROJECT_ROOT, args.data006_config)
        if not os.path.isfile(cfg):
            sys.exit("ERROR: DATA-006 config not found: %s" % cfg)
    if not os.path.isfile(os.path.join(PROJECT_ROOT, args.mainconf)):
        sys.exit("ERROR: mainconf not found: %s" % args.mainconf)
    if args.three_event:
        if args.mode != "phi":
            sys.exit("ERROR: --three-event is only defined for --mode phi")
        if not args.three_event_cap or args.three_event_cap <= 0:
            sys.exit("ERROR: --three-event requires --three-event-cap > 0")
        if not args.three_event_seed or args.three_event_seed <= 0:
            sys.exit("ERROR: --three-event requires --three-event-seed > 0")

    blocks = read_blocks(args)
    counts = (block_event_counts(blocks, os.path.join(workdir, "block_events.txt"))
              if args.events_per_job > 0 else {})
    build_runtime(bundle, args.reuse_runtime)

    units = []
    for tree in blocks:
        if args.events_per_job > 0:
            total = counts[tree]
            nslice = (total + args.events_per_job - 1) // args.events_per_job
            for islice in range(nslice):
                units.append((tree, "%s_s%02d" % (block_key(tree), islice),
                              islice * args.events_per_job, args.events_per_job))
        else:
            units.append((tree, block_key(tree), 0, args.max_events))

    submitted, skipped, condor = [], [], []
    for tree, key, first_event, max_events in units:
        name = "downstream_%s_%s" % (args.tag, key)
        out_root = os.path.join(outdir, "%s.root" % name)
        if not args.force and os.path.isfile(out_root) and os.path.getsize(out_root) > 0:
            skipped.append(name)
            continue
        if args.limit and len(submitted) >= args.limit:
            continue
        sh = os.path.join(workdir, "job", "%s.sh" % name)
        with open(sh, "w") as fh:
            fh.write(job_script(args, bundle, tree, out_root, first_event, max_events))
        os.chmod(sh, 0o755)
        condor.append(condor_block(args, workdir, name))
        submitted.append(name)

    if not submitted:
        print("[submit] nothing to do: %d block(s) already have output" % len(skipped))
        return 0

    submit_file = os.path.join(workdir, "job", "downstream_%s.condor" % args.tag)
    with open(submit_file, "w") as fh:
        fh.write("".join(condor))
    print("[submit] mode=%s blocks=%d jobs=%d submit=%d skip=%d" %
          (args.mode, len(blocks), len(units), len(submitted), len(skipped)))
    if args.three_event:
        print("[submit] three-event cap=%d seed=%d" %
              (args.three_event_cap, args.three_event_seed))
    print("[submit] outputs  -> %s" % outdir)
    print("[submit] submit file -> %s" % submit_file)
    if args.dry_run:
        print("[submit] --dry-run: not calling condor_submit")
        return 0
    rc = subprocess.call(["condor_submit", submit_file])
    if rc != 0:
        sys.exit("ERROR: condor_submit exited %d" % rc)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
