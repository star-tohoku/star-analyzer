#!/usr/bin/env python3
"""Build the Step 7 completeness ledger for the reduced-tree production.

The plan's Step 7 requires that every input PicoDst end up either "processed" or "excluded with a
reason", with no unaccounted input, no duplicates and no ambiguous provenance. That means joining
five records that nobody joins for you:

  1. the SUMS session   -- which process was given which PicoDst, and how many events the catalog
                           says each one holds (the per-process .list files in the artifacts tar)
  2. the job stderr     -- what happened to a process that left no output
  3. the merged blocks  -- which subjob outputs were merged (the .inputs list beside each block)
  4. SourceFileTable    -- which PicoDst a surviving output actually read
  5. hCounters          -- how many events those reads delivered

Items 4 and 5 come from tools/treecheck/SourceLedger.C, which dumps them per block; this script
does the joining and the classification. It only reads, and writes TSV.

Usage
  script/step7_ledger.py --sums <dir with the unpacked SUMS artifacts> \
                         --err <dir with stderr.<requestId>_<proc>.err> \
                         --merged <dir with *.inputs beside every merged block> \
                         --src <SourceLedger .src.tsv> --counters <SourceLedger .cnt.tsv> \
                         --out <output prefix>

Writes <out>.tsv (one row per input PicoDst) and prints the summary tables.
"""

import argparse
import glob
import os
import re
import sys
from collections import defaultdict

PICO = re.compile(r"st_physics_\d+_raw_\d+\.picoDst\.root")


def read_sums_lists(sums_dir):
    """proc -> [(fileName, url, catalogEvents)], from the per-process .list files."""
    per_proc = {}
    pattern = os.path.join(sums_dir, "*_[0-9]*.list")
    for path in glob.glob(pattern):
        m = re.search(r"_(\d+)\.list$", os.path.basename(path))
        if not m:
            continue
        proc = int(m.group(1))
        rows = []
        for line in open(path):
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            url = parts[0]
            events = int(parts[1]) if len(parts) > 1 else -1
            rows.append((os.path.basename(url), url, events))
        per_proc[proc] = rows
    if not per_proc:
        sys.exit("ERROR: no per-process .list files under %s" % sums_dir)
    return per_proc


def read_merged_inputs(merged_dir):
    """The set of process indices whose output reached a merged block."""
    produced = set()
    blocks = {}
    for path in sorted(glob.glob(os.path.join(merged_dir, "*.inputs"))):
        block = os.path.basename(path)[:-len(".inputs")] + ".root"
        for line in open(path):
            line = line.strip()
            if not line:
                continue
            m = re.search(r"_(\d+)\.root$", line)
            if not m:
                sys.exit("ERROR: cannot read a process index from %s" % line)
            proc = int(m.group(1))
            if proc in produced:
                sys.exit("ERROR: process %d appears in more than one merged block" % proc)
            produced.add(proc)
            blocks[proc] = block
    if not produced:
        sys.exit("ERROR: no *.inputs with content under %s" % merged_dir)
    return produced, blocks


def read_source_table(src_tsv):
    """(proc, fileName) -> block, from the merged SourceFileTable rows.

    subjobId is the SUMS process index plus one (plan sec 13.1), which is what makes this join
    possible at all: sourceFileIndex alone is local to one output file.
    """
    recorded = {}
    with open(src_tsv) as fh:
        header = fh.readline()
        if not header.startswith("block\t"):
            sys.exit("ERROR: %s does not look like a SourceLedger .src.tsv" % src_tsv)
        for line in fh:
            block, subjob, _index, name = line.rstrip("\n").split("\t")
            key = (int(subjob) - 1, name)
            if key in recorded:
                sys.exit("ERROR: %s read twice by process %d" % (name, int(subjob) - 1))
            recorded[key] = block
    return recorded


def classify_no_output(err_dir, request_id, proc):
    """Why a process left no output, from its stderr alone."""
    path = os.path.join(err_dir, "stderr.%s_%d.err" % (request_id, proc))
    if not os.path.isfile(path):
        return "no_stderr"
    if os.path.getsize(path) == 0:
        return "never_started"
    text = open(path, errors="replace").read()
    if "cannot stat" in text:
        return "input_copy_failed"
    if "Terminated" in text:
        return "terminated"
    return "unknown"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sums", required=True)
    ap.add_argument("--err", required=True)
    ap.add_argument("--merged", required=True)
    ap.add_argument("--src", required=True)
    ap.add_argument("--counters", required=True)
    ap.add_argument("--request-id", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    for path in (args.sums, args.err, args.merged):
        if not os.path.isdir(path):
            sys.exit("ERROR: not a directory: %s" % path)
    for path in (args.src, args.counters):
        if not os.path.isfile(path):
            sys.exit("ERROR: not a file: %s" % path)

    per_proc = read_sums_lists(args.sums)
    produced, block_of = read_merged_inputs(args.merged)
    recorded = read_source_table(args.src)

    # A job that read nothing at all still writes an output, an empty SourceFileTable and a zero
    # hCounters: the chain macro calls Finish() on "no entries found". Those outputs are real but
    # carry no input, and the ledger has to say so rather than count them as processed.
    read_by_proc = defaultdict(int)
    for (proc, _name) in recorded:
        read_by_proc[proc] += 1

    rows = []
    counts = defaultdict(int)
    events = defaultdict(int)
    for proc in sorted(per_proc):
        no_output_reason = None
        if proc not in produced:
            no_output_reason = classify_no_output(args.err, args.request_id, proc)
        for name, url, nev in per_proc[proc]:
            scheme = "xrootd" if url.startswith("root:") else "nfs_path"
            block = recorded.get((proc, name))
            if block:
                status, reason = "processed", ""
            elif no_output_reason:
                status, reason = "excluded_no_output", no_output_reason
            elif read_by_proc[proc] == 0:
                status, reason = "excluded_unread", "job_read_no_input"
            else:
                status, reason = "excluded_unread", "staging_failed_" + scheme
            rows.append((name, proc, proc + 1, scheme, nev, status, reason, block or ""))
            counts[(status, reason)] += 1
            events[(status, reason)] += nev

    out_tsv = args.out + ".tsv"
    with open(out_tsv, "w") as fh:
        fh.write("fileName\tprocId\tsubjobId\tscheme\tcatalogEvents\tstatus\treason\tblock\n")
        for row in rows:
            fh.write("\t".join(str(x) for x in row) + "\n")

    total_files = len(rows)
    processed = sum(n for (s, _r), n in counts.items() if s == "processed")
    print("input PicoDst      %d" % total_files)
    print("processed          %d" % processed)
    print("excluded           %d" % (total_files - processed))
    print()
    print("%-22s %-24s %8s %16s" % ("status", "reason", "files", "catalogEvents"))
    for (status, reason), n in sorted(counts.items()):
        print("%-22s %-24s %8d %16d" % (status, reason, n, events[(status, reason)]))

    # The completeness criterion: every input is processed or excluded with a reason, nothing is
    # counted twice, and the ledger's processed count is the number of SourceFileTable rows.
    print()
    print("unaccounted        %d" % (total_files - sum(counts.values())))
    print("SourceFileTable    %d rows  (ledger processed: %d)" % (len(recorded), processed))
    print("duplicate reads    %d" % (len(recorded) - len(set(k[1] for k in recorded))))

    delivered = 0
    with open(args.counters) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = header.index("nInputEvents")
        for line in fh:
            delivered += int(line.rstrip("\n").split("\t")[col])
    expected = sum(nev for (_n, _p, _s, _sc, nev, st, _r, _b) in
                   [(r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7]) for r in rows] if st == "processed")
    print("catalog events on processed input  %d" % expected)
    print("hCounters nInputEvents             %d" % delivered)
    print("difference                         %d (%.2f %%)"
          % (delivered - expected, 100.0 * (delivered - expected) / expected))
    print()
    print("ledger written to %s" % out_tsv)
    return 0


if __name__ == "__main__":
    sys.exit(main())
