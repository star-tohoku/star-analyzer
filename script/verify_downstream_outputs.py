#!/usr/bin/env python3
"""Check that a `submit_downstream_blocks.py` run really finished.

The submitter skips a slice whose output file exists and is not empty, and ROOT creates that file
the moment a job starts, so file existence proves nothing. This checks the job logs instead:

* every `log/*.out` ends in `exit 0`;
* every job's output ROOT file exists and is not empty;
* the processed-event counts add up to the block totals in `block_events.txt`, when the run was
  sliced (`--events-per-job`).

It prints what is missing and, with `--print-resubmit`, the job names to resubmit with
`--force`.
"""

import argparse
import glob
import os
import re
import sys

PROJECT_ROOT = os.path.realpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
EXIT_RE = re.compile(r"^exit (\d+)\b", re.M)
EVENTS_RE = re.compile(r"\[downstreamV3\] events=(\d+)/(\d+) \[(\d+),(\d+)\)")
EVENTS_RE_OLD = re.compile(r"\[downstreamV3\] events=(\d+)/(\d+)")
SLICE_SUFFIX_RE = re.compile(r"_s\d+$")


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--workdir", default=None)
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--print-resubmit", action="store_true")
    args = ap.parse_args(argv)

    workdir = os.path.realpath(args.workdir or os.path.join(PROJECT_ROOT, "job/merge", args.tag))
    outdir = os.path.realpath(
        args.outdir or "/gpfs/mnt/gpfs01/star/pwg/oura/rootfile/%s" % args.tag)

    jobs = sorted(glob.glob(os.path.join(workdir, "job", "downstream_*.sh")))
    if not jobs:
        sys.exit("ERROR: no job scripts under %s/job" % workdir)

    ok, failed, unfinished, empty = [], [], [], []
    events_by_block = {}
    for job in jobs:
        name = os.path.basename(job)[:-len(".sh")]
        out_log = os.path.join(workdir, "log", "%s.out" % name)
        root_out = os.path.join(outdir, "%s.root" % name)
        if not os.path.isfile(out_log):
            unfinished.append(name)
            continue
        with open(out_log, errors="replace") as fh:
            text = fh.read()
        codes = EXIT_RE.findall(text)
        if not codes:
            unfinished.append(name)
            continue
        if codes[-1] != "0":
            failed.append("%s (exit %s)" % (name, codes[-1]))
            continue
        if not os.path.isfile(root_out) or os.path.getsize(root_out) == 0:
            empty.append(name)
            continue
        m = EVENTS_RE.search(text) or EVENTS_RE_OLD.search(text)
        if m:
            block = SLICE_SUFFIX_RE.sub("", name)
            events_by_block.setdefault(block, [0, int(m.group(2))])[0] += int(m.group(1))
        ok.append(name)

    print("[verify] tag=%s jobs=%d ok=%d failed=%d unfinished=%d empty-output=%d"
          % (args.tag, len(jobs), len(ok), len(failed), len(unfinished), len(empty)))
    for label, items in (("FAILED", failed), ("UNFINISHED", unfinished), ("EMPTY", empty)):
        for item in items[:20]:
            print("  %-11s %s" % (label, item))
        if len(items) > 20:
            print("  %-11s ... and %d more" % (label, len(items) - 20))

    short = [(b, got, total) for b, (got, total) in sorted(events_by_block.items())
             if got != total]
    if short:
        print("[verify] %d block(s) do not add up to their event total:" % len(short))
        for b, got, total in short[:20]:
            print("  SHORT       %s %d/%d" % (b, got, total))

    if args.print_resubmit and (failed or unfinished or empty):
        print("[verify] resubmit these outputs (delete them, then rerun the submitter):")
        for name in failed + unfinished + empty:
            print(os.path.join(outdir, "%s.root" % name.split(" ")[0]))

    return 0 if not (failed or unfinished or empty or short) else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
