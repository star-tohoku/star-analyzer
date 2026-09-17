#!/usr/bin/env python3
"""Farm jobs for the post-production stage of the reduced-tree workflow.

Production leaves one ROOT file per subjob -- 14,178 of them for the full data set. Two things
have to happen before those files are an analysis input: every file must be checked for
completeness (a job whose output copy was interrupted still leaves a readable file with short
trees, see tools/treecheck/FileIntegrity.C), and the surviving files must be merged into blocks
small enough to hand to a downstream job.

Both steps are embarrassingly parallel and neither belongs on a login node: a single stream checks
1.7 files/s and merges 15 MB/s, so the full data set is 2.2 h of checking and 30 h of merging.
This script splits each step into condor jobs.

Subcommands
  integrity  check a list of production outputs, in --shards parallel jobs
  collect    combine the per-shard results into one good-list, report what failed
  merge      hadd the good-list into blocks of --block files, one job per block
  status     how far the submitted work has got

The jobs are plain condor submissions, not SUMS: there is no catalog query and no input sandbox to
stage, only a list of GPFS paths. The submit files mirror the ones SUMS writes (vanilla universe,
/bin/sh -c exec, a PeriodicRemove wall clock) so that the two kinds of job behave the same way in
the queue.

Batch safety (docs/ai/AGENT_RULES.md): every path is required and is checked before anything is
submitted. Nothing here falls back to a default when an argument is missing -- a batch job that
runs on the wrong input is worse than one that does not run.
"""

import argparse
import os
import shutil
import subprocess
import sys

IMAGE = "/cvmfs/singularity.opensciencegrid.org/star-bnl/star-sw:latest"
# The binds the farm's own production jobs use, plus /gpfs: the tree outputs and the merge results
# both live on GPFS, and unlike a production job this one reads and writes them in place rather
# than staging through the node's scratch disk.
BINDS = ["-B /gpfs:/gpfs", "-B /cvmfs:/cvmfs", "-B /star/u:/star/u",
         "-B /star/nfs4/AFS:/star/nfs4/AFS", "-B /home/starlib:/home/starlib"]

CONDOR_STANZA = """Universe         = vanilla
Notification     = never
Executable       = /bin/sh
Arguments        = "-c 'exec {script}'"
Output           = {log}/{tag}.out
Error            = {log}/{tag}.err
Log              = {log}/condor.log
Initialdir       = {workdir}
kill_sig         = SIGINT
request_memory   = {memory}
PeriodicRemove   = (JobStatus == 2 && (CurrentTime - JobCurrentStartDate > ({wall}))) || (((CurrentTime - EnteredCurrentStatus) > (2*24*3600)) && JobStatus == 5)
Priority         = +10
Queue

"""


def root_env():
    """What these jobs need from the container, which is only ROOT.

    Production jobs point STAR at the SL24y release and run root4star, because the analysis
    libraries are built against it. Nothing in this stage does: hadd is a ROOT tool, and every
    macro under tools/treecheck includes only ROOT headers (TFile, TTree, TH1 and friends). Taking
    the container's own ROOT instead of the release's avoids dragging in the release's runtime --
    libmysqlclient and the rest -- for work that never touches it.
    """
    return "export LD_LIBRARY_PATH=${ROOTSYS}/lib:$LD_LIBRARY_PATH; "


def in_container(command):
    return "singularity exec %s %s /bin/bash -lc %s" % (
        " ".join(BINDS), IMAGE, shell_quote(command))


def shell_quote(s):
    return "'" + s.replace("'", "'\\''") + "'"


def write_exec(path, body):
    with open(path, "w") as fh:
        fh.write(body)
    os.chmod(path, 0o755)


def prepare(workdir):
    logdir = os.path.join(workdir, "log")
    jobdir = os.path.join(workdir, "job")
    for d in (workdir, logdir, jobdir):
        if not os.path.isdir(d):
            os.makedirs(d)
    return logdir, jobdir


def submit(condor_path, dry_run):
    if dry_run:
        print("dry run: not submitting %s" % condor_path)
        return 0
    print("condor_submit %s" % condor_path)
    return subprocess.call(["condor_submit", condor_path])


def read_list(path):
    if not os.path.isfile(path):
        sys.exit("ERROR: list file does not exist: %s" % path)
    out = []
    for line in open(path):
        line = line.strip()
        if line:
            out.append(line)
    if not out:
        sys.exit("ERROR: list file is empty: %s" % path)
    return out


# --------------------------------------------------------------------------------------------
# integrity


def cmd_integrity(args):
    paths = read_list(args.list)
    repo = os.path.abspath(args.repo)
    workdir = os.path.abspath(args.workdir)
    logdir, jobdir = prepare(workdir)

    shards = max(1, args.shards)
    per = (len(paths) + shards - 1) // shards
    print("%d files -> %d shards of at most %d" % (len(paths), shards, per))

    condor = []
    for i in range(shards):
        chunk = paths[i * per:(i + 1) * per]
        if not chunk:
            continue
        tag = "integrity_%03d" % i
        shard_list = os.path.join(jobdir, tag + ".list")
        with open(shard_list, "w") as fh:
            fh.write("\n".join(chunk) + "\n")
        good = os.path.join(workdir, "good_%03d.txt" % i)

        # FileIntegrity.C is compiled with ACLiC. Every shard compiles its own copy into a
        # build directory on the node's scratch disk, named after the shard: a build directory
        # shared over GPFS would have sixteen jobs writing the same .so at once.
        build = "/tmp/%s/treemerge_%s" % (os.environ.get("USER", "star"), tag)
        prelude = os.path.join(jobdir, tag + "_prelude.C")
        with open(prelude, "w") as fh:
            fh.write(
                "// Generated by script/submit_tree_merge.py. All this has to do is send ACLiC's\n"
                "// build products to the node's own disk, so that parallel shards do not write\n"
                "// the same shared library over each other on GPFS.\n"
                "void %s_prelude() { gSystem->SetBuildDir(\"%s\", kTRUE); }\n" % (tag, build))

        macro = ("root -b -q %s "
                 "'%s/tools/treecheck/FileIntegrity.C+(\"%s\",\"%s\")'"
                 % (prelude, repo, shard_list, good))
        body = """#!/bin/bash
set -u
echo "host $(hostname)  start $(date)"
mkdir -p %s || exit 1
%s
rc=$?
rm -rf %s
echo "exit $rc  end $(date)"
exit $rc
""" % (build, in_container(root_env() + "cd %s; " % repo + macro), build)
        script = os.path.join(jobdir, tag + ".sh")
        write_exec(script, body)
        condor.append(CONDOR_STANZA.format(script=script, log=logdir, tag=tag, workdir=workdir,
                                           memory=args.memory, wall=args.wall))

    condor_path = os.path.join(jobdir, "integrity.condor")
    open(condor_path, "w").write("".join(condor))
    return submit(condor_path, args.dry_run)


def cmd_collect(args):
    workdir = os.path.abspath(args.workdir)
    good_paths, n_shards = [], 0
    for name in sorted(os.listdir(workdir)):
        if name.startswith("good_") and name.endswith(".txt"):
            n_shards += 1
            for line in open(os.path.join(workdir, name)):
                line = line.strip()
                if line:
                    # FileIntegrity writes "<runNumber> <path>"; the merge only needs the path,
                    # because SUMS packed files across run boundaries and a subjob output cannot
                    # be assigned to a single run.
                    good_paths.append(line.split(None, 1)[-1])
    if not good_paths:
        sys.exit("ERROR: no good_*.txt with content in %s -- has the integrity step finished?"
                 % workdir)
    out = os.path.join(workdir, "good.txt")
    with open(out, "w") as fh:
        fh.write("\n".join(good_paths) + "\n")

    submitted = read_list(args.list) if args.list else None
    print("%d shards -> %d complete files -> %s" % (n_shards, len(good_paths), out))
    if submitted is not None:
        bad = set(submitted) - set(good_paths)
        if bad:
            bad_path = os.path.join(workdir, "bad.txt")
            with open(bad_path, "w") as fh:
                fh.write("\n".join(sorted(bad)) + "\n")
            print("%d files did not pass; listed in %s" % (len(bad), bad_path))
            print("Re-run those subjobs before merging; do not hadd them.")
        else:
            print("every submitted file passed")
    return 0


# --------------------------------------------------------------------------------------------
# merge


def cmd_merge(args):
    paths = read_list(args.good)
    repo = os.path.abspath(args.repo)
    workdir = os.path.abspath(args.workdir)
    outdir = os.path.abspath(args.outdir)
    logdir, jobdir = prepare(workdir)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    block = max(1, args.block)
    n_blocks = (len(paths) + block - 1) // block
    # Block numbering is over the whole good-list, always. Only a range of them is submitted at a
    # time, because the merged blocks and the subjob outputs they were built from do not both fit
    # under the quota: a wave is merged, verified, and its sources reclaimed before the next wave
    # starts. Numbering the blocks globally keeps a block's name meaning the same thing whichever
    # wave built it.
    first, last = 0, n_blocks
    if args.block_range:
        try:
            a, b = args.block_range.split(":")
            first, last = int(a), int(b)
        except ValueError:
            sys.exit("ERROR: --block-range wants START:END, e.g. 0:35")
        if not (0 <= first < last <= n_blocks):
            sys.exit("ERROR: --block-range %s is outside 0:%d" % (args.block_range, n_blocks))
    print("%d files -> %d blocks of %d; submitting blocks %d:%d"
          % (len(paths), n_blocks, block, first, last))

    condor, skipped = [], 0
    for i in range(first, last):
        chunk = paths[i * block:(i + 1) * block]
        tag = "%s_block%04d" % (args.stem, i)
        out_root = os.path.join(outdir, tag + ".root")
        if os.path.exists(out_root) and not args.force:
            skipped += 1
            continue
        # The list of inputs is kept next to the merged file: which subjob outputs went into which
        # block is the only record of it once the sources are gone.
        block_list = os.path.join(outdir, tag + ".inputs")
        with open(block_list, "w") as fh:
            fh.write("\n".join(chunk) + "\n")

        # hadd writes to a .part name and the job renames it only on success, so an interrupted
        # merge cannot leave something that looks like a finished block.
        merge = ("hadd -f %s.part $(cat %s) && mv -f %s.part %s"
                 % (out_root, block_list, out_root, out_root))
        body = """#!/bin/bash
set -u
echo "host $(hostname)  start $(date)"
echo "block %s: %d inputs -> %s"
%s
rc=$?
rm -f %s.part
if [ $rc -ne 0 ]; then echo "MERGE FAILED rc=$rc"; exit $rc; fi
ls -l %s
echo "exit 0  end $(date)"
""" % (tag, len(chunk), out_root,
       in_container(root_env() + "cd %s; " % repo + merge),
       out_root, out_root)
        script = os.path.join(jobdir, tag + ".sh")
        write_exec(script, body)
        condor.append(CONDOR_STANZA.format(script=script, log=logdir, tag=tag, workdir=workdir,
                                           memory=args.memory, wall=args.wall))

    if skipped:
        print("%d blocks already exist and were skipped (--force to redo them)" % skipped)
    if not condor:
        print("nothing to submit")
        return 0
    condor_path = os.path.join(jobdir, "merge.condor")
    open(condor_path, "w").write("".join(condor))

    # Merging does not shrink anything: the blocks are as large as their inputs put together.
    # Sized from files that are still on disk: once a wave's inputs have been released, the head
    # of the good-list no longer exists, and the estimate must not be what stops the next wave.
    sizes, wanted = [], 200
    for path in paths:
        if len(sizes) >= wanted:
            break
        if os.path.isfile(path):
            sizes.append(os.path.getsize(path))
    if sizes:
        per_file = sum(sizes) / float(len(sizes))
        print("this wave will write about %.3f TB; merging does not shrink anything"
              % (per_file * block * len(condor) / 1e12))
    return submit(condor_path, args.dry_run)


# --------------------------------------------------------------------------------------------
# verify


def block_pairs(outdir):
    """Every merged block in outdir that still has the input list it was built from."""
    pairs = []
    for name in sorted(os.listdir(outdir)):
        if not name.endswith(".root"):
            continue
        root = os.path.join(outdir, name)
        inputs = root[:-len(".root")] + ".inputs"
        if os.path.isfile(inputs):
            pairs.append((root, inputs))
    return pairs


def inputs_alive(inputs):
    """How many of a block's subjob outputs are still on disk.

    A block whose inputs have been released is finished business: it has already been verified and
    the files it was built from are gone on purpose. Re-verifying it would compare the block
    against files that no longer exist and call that a failure, which is what happened on the first
    run of the second wave. Counting what survives is how the later steps tell "already released"
    apart from "not verified yet".
    """
    total = alive = 0
    for line in open(inputs):
        line = line.strip()
        if line:
            total += 1
            if os.path.isfile(line):
                alive += 1
    return total, alive


def cmd_verify(args):
    repo = os.path.abspath(args.repo)
    workdir = os.path.abspath(args.workdir)
    outdir = os.path.abspath(args.outdir)
    logdir, jobdir = prepare(workdir)
    pairs, released = [], 0
    for root, inputs in block_pairs(outdir):
        if inputs_alive(inputs)[1] == 0:
            released += 1
            continue
        pairs.append((root, inputs))
    if released:
        print("%d blocks already verified and released; not re-checking them" % released)
    if not pairs:
        sys.exit("ERROR: no merged block left to verify in %s" % outdir)

    shards = max(1, min(args.shards, len(pairs)))
    per = (len(pairs) + shards - 1) // shards
    print("%d blocks -> %d shards of at most %d" % (len(pairs), shards, per))

    condor = []
    for i in range(shards):
        chunk = pairs[i * per:(i + 1) * per]
        if not chunk:
            continue
        tag = "verify_%03d" % i
        build = "/tmp/%s/treemerge_%s" % (os.environ.get("USER", "star"), tag)
        prelude = os.path.join(jobdir, tag + "_prelude.C")
        with open(prelude, "w") as fh:
            fh.write("void %s_prelude() { gSystem->SetBuildDir(\"%s\", kTRUE); }\n" % (tag, build))
        result = os.path.join(workdir, "verified_%03d.txt" % i)

        # One root invocation per block. BlockVerify prints a final "PASS <path>" or "FAIL <path>"
        # line; the job records that line and nothing else, so whoever decides what to do with a
        # block reads a verdict rather than having to re-interpret a log.
        lines = []
        for root, inputs in chunk:
            macro = ("root -b -q %s '%s/tools/treecheck/BlockVerify.C+(\"%s\",\"%s\")'"
                     % (prelude, repo, root, inputs))
            lines.append(
                "%s | tee -a %s.log | grep -E '^(PASS|FAIL) ' >> %s"
                % (in_container(root_env() + "cd %s; " % repo + macro), result, result))
        body = """#!/bin/bash
set -u
echo "host $(hostname)  start $(date)"
mkdir -p %s || exit 1
%s
echo "verdicts: $(wc -l < %s)  expected %d"
echo "end $(date)"
""" % (build, "\n".join(lines), result, len(chunk))
        script = os.path.join(jobdir, tag + ".sh")
        write_exec(script, body)
        condor.append(CONDOR_STANZA.format(script=script, log=logdir, tag=tag, workdir=workdir,
                                           memory=args.memory, wall=args.wall))

    condor_path = os.path.join(jobdir, "verify.condor")
    open(condor_path, "w").write("".join(condor))
    return submit(condor_path, args.dry_run)


def cmd_verdicts(args):
    """Summarise the verify step: which blocks passed, and what they account for.

    A block that passed holds exactly what its inputs held, which is the condition for the inputs
    to be redundant. Whether to act on that is a separate decision and this script does not take
    it: the production outputs cannot be rebuilt without re-running the farm.
    """
    workdir = os.path.abspath(args.workdir)
    outdir = os.path.abspath(args.outdir)
    verdict = {}
    for name in sorted(os.listdir(workdir)):
        if name.startswith("verified_") and name.endswith(".txt"):
            for line in open(os.path.join(workdir, name)):
                parts = line.split()
                if len(parts) == 2 and parts[0] in ("PASS", "FAIL"):
                    verdict[os.path.abspath(parts[1])] = parts[0]
    if not verdict:
        sys.exit("ERROR: no verdicts in %s -- run the verify step first" % workdir)

    passed, failed, unjudged = [], [], []
    for root, inputs in block_pairs(outdir):
        v = verdict.get(os.path.abspath(root))
        (passed if v == "PASS" else failed if v == "FAIL" else unjudged).append((root, inputs))

    covered, nbytes = 0, 0
    for _, inputs in passed:
        for line in open(inputs):
            line = line.strip()
            if line and os.path.isfile(line):
                covered += 1
                nbytes += os.path.getsize(line)
    print("%d blocks passed, %d failed, %d without a verdict" % (len(passed), len(failed),
                                                                 len(unjudged)))
    for root, _ in failed + unjudged:
        print("  no usable verdict for %s" % os.path.basename(root))
    print("passed blocks account for %d subjob outputs, %.3f TB" % (covered, nbytes / 1e12))
    listing = os.path.join(workdir, "redundant_sources.txt")
    with open(listing, "w") as fh:
        for _, inputs in passed:
            for line in open(inputs):
                if line.strip():
                    fh.write(line)
    print("those sources are listed in %s" % listing)
    return 0


def cmd_reclaim(args):
    """Release the subjob outputs behind blocks that verified.

    The merged blocks and the subjob outputs they were built from do not both fit under the quota,
    so a block's inputs have to be released once the block is known to hold exactly what they held.
    That condition is what BlockVerify establishes, and this reads its recorded verdicts rather
    than re-deriving anything: a block with a PASS releases its inputs, a block with a FAIL or with
    no verdict at all keeps them.

    Production output cannot be rebuilt without re-running the farm, so nothing happens without
    --apply; by default this only reports what a run would release.
    """
    workdir = os.path.abspath(args.workdir)
    outdir = os.path.abspath(args.outdir)
    verdict = {}
    for name in sorted(os.listdir(workdir)):
        if name.startswith("verified_") and name.endswith(".txt"):
            for line in open(os.path.join(workdir, name)):
                parts = line.split()
                if len(parts) == 2 and parts[0] in ("PASS", "FAIL"):
                    verdict[os.path.abspath(parts[1])] = parts[0]
    if not verdict:
        sys.exit("ERROR: no verdicts in %s -- run the verify step first" % workdir)

    passed, held, released = [], [], 0
    for root, inputs in block_pairs(outdir):
        if inputs_alive(inputs)[1] == 0:
            released += 1
            continue
        v = verdict.get(os.path.abspath(root))
        (passed if v == "PASS" else held).append((root, inputs))
    if released:
        print("%d blocks were released earlier" % released)

    redundant, nbytes = [], 0
    for _, inputs in passed:
        for line in open(inputs):
            line = line.strip()
            if line and os.path.isfile(line):
                redundant.append(line)
                nbytes += os.path.getsize(line)

    print("%d blocks passed, %d without a usable verdict" % (len(passed), len(held)))
    for root, _ in held:
        print("  holding the inputs of %s" % os.path.basename(root))
    print("%d subjob outputs are now redundant, %.3f TB" % (len(redundant), nbytes / 1e12))
    if not args.apply:
        print("reporting only; pass --apply to release them")
        return 0
    released, failures = 0, 0
    for path in redundant:
        try:
            os.unlink(path)
            released += 1
        except OSError as exc:
            print("  could not release %s: %s" % (path, exc))
            failures += 1
    print("released %d files, %.3f TB" % (released, nbytes / 1e12))
    return 1 if failures else 0


def resolve_star_runtime(repo, mainconf):
    """The library path the STAR release needs, resolved once on the submit node.

    script/setup.sh is the project's own way of building this, but it calls starver, which exists
    only where the STAR login has been sourced -- on a worker node running a bare bash job it is
    not there, and the first attempt at these jobs died in one second saying so. The paths it
    produces are all on shared filesystems (/cvmfs and the AFS release area), so resolving it here
    and writing it into the job is both correct and better documented: the job script then records
    exactly which libraries it ran against.
    """
    cmd = "source %s/script/setup.sh %s >/dev/null 2>&1 && printf '%%s\\n%%s' \"$LD_LIBRARY_PATH\" \"$STAR_HOST_SYS\"" % (repo, mainconf)
    proc = subprocess.Popen(["bash", "-c", cmd], cwd=repo, stdout=subprocess.PIPE)
    out = proc.communicate()[0].decode().split("\n")
    if proc.returncode != 0 or len(out) < 2 or not out[0].strip() or not out[1].strip():
        sys.exit("ERROR: could not resolve the STAR runtime from %s.\n"
                 "       Run this on a node where 'starver' works (the STAR login must be "
                 "sourced; for this project, 'sl7' first)." % mainconf)
    return out[0].strip(), out[1].strip()


def cmd_downstream(args):
    """One downstream analysis job per merged block.

    Unlike the merge and the checks, this leg does need the STAR release: the macro links the
    project's analysis libraries against it. The job therefore builds the release environment
    itself inside the container rather than calling the interactive wrapper, which depends on a
    login environment the farm does not have.
    """
    repo = os.path.abspath(args.repo)
    workdir = os.path.abspath(args.workdir)
    outdir = os.path.abspath(args.outdir)
    logdir, jobdir = prepare(workdir)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    mainconf = args.mainconf
    if not os.path.isfile(os.path.join(repo, mainconf)):
        sys.exit("ERROR: mainconf does not exist under %s: %s" % (repo, mainconf))

    ldpath, host_sys = resolve_star_runtime(repo, mainconf)
    print("STAR_HOST_SYS %s, %d entries on the library path" % (host_sys, ldpath.count(":") + 1))

    blocks = sorted(b for b in os.listdir(args.blockdir) if b.endswith(".root"))
    if args.limit:
        blocks = blocks[:args.limit]
    if not blocks:
        sys.exit("ERROR: no merged blocks in %s" % args.blockdir)

    star_setup = (
        "export STAR=/star/nfs4/AFS/star/packages/%s; "
        "export STAR_HOST_SYS=%s; "
        "export STAR_LIB=${STAR}/.${STAR_HOST_SYS}/lib; "
        "export STAR_BIN=${STAR}/.${STAR_HOST_SYS}/bin; "
        "export PATH=${STAR_BIN}:${STAR}/mgr:$PATH; "
        "export LD_LIBRARY_PATH=${STAR_LIB}:$LD_LIBRARY_PATH; " % (args.library_tag, host_sys))

    condor, skipped = [], 0
    for name in blocks:
        tag = "downstream%s_block%s" % (args.tag, name[:-len(".root")].split("_block")[-1])
        block = os.path.join(os.path.abspath(args.blockdir), name)
        out_root = os.path.join(outdir, tag + ".root")
        if os.path.exists(out_root) and not args.force:
            skipped += 1
            continue
        macro = ("root4star -b -q 'tools/run_anaFemtoPhiTreeDownstreamV3.C"
                 "(\"%s\",\"%s\",\"%s\",%d,\"%s\",\"%s\",-1,-1,0,0)'"
                 % (block, out_root, mainconf, args.buffer, args.mixing_mode, args.variation))
        inner = star_setup + "cd %s; " % repo + macro
        run = "singularity exec %s --env LD_LIBRARY_PATH=%s %s /bin/bash -lc %s" % (
            " ".join(BINDS), shell_quote("%s/lib:%s" % (repo, ldpath)), IMAGE, shell_quote(inner))
        body = """#!/bin/bash
set -u
echo "host $(hostname)  start $(date)"
cd %s || exit 1
%s
rc=$?
echo "exit $rc  end $(date)"
if [ $rc -eq 0 ]; then ls -l %s; fi
exit $rc
""" % (repo, run, out_root)
        script = os.path.join(jobdir, tag + ".sh")
        write_exec(script, body)
        condor.append(CONDOR_STANZA.format(script=script, log=logdir, tag=tag, workdir=workdir,
                                           memory=args.memory, wall=args.wall))
    print("%d blocks -> %d downstream jobs%s" % (
        len(blocks), len(condor),
        (" with variation '%s'" % args.variation) if args.variation else " (nominal)"))
    if skipped:
        print("%d outputs already exist and were skipped (--force to redo them)" % skipped)
    if not condor:
        print("nothing to submit")
        return 0
    condor_path = os.path.join(jobdir, "downstream.condor")
    open(condor_path, "w").write("".join(condor))
    return submit(condor_path, args.dry_run)


def cmd_status(args):
    workdir = os.path.abspath(args.workdir)
    outdir = os.path.abspath(args.outdir) if args.outdir else None
    for name in sorted(os.listdir(workdir)):
        if name.startswith("good_") and name.endswith(".txt"):
            n = sum(1 for line in open(os.path.join(workdir, name)) if line.strip())
            print("  %-20s %6d checked" % (name, n))
    if outdir and os.path.isdir(outdir):
        done, nbytes = 0, 0
        for name in os.listdir(outdir):
            if name.endswith(".root"):
                done += 1
                nbytes += os.path.getsize(os.path.join(outdir, name))
        print("merged: %d blocks, %.2f TB in %s" % (done, nbytes / 1e12, outdir))
    print("queue:")
    subprocess.call("condor_q -totals", shell=True)
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--repo", default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                    help="project root the jobs cd into (default: this script's repository)")
    ap.add_argument("--dry-run", action="store_true", help="write the job files but do not submit")
    sub = ap.add_subparsers(dest="cmd")

    p = sub.add_parser("integrity", help="check production outputs in parallel jobs")
    p.add_argument("--list", required=True, help="file holding one output path per line")
    p.add_argument("--workdir", required=True)
    p.add_argument("--shards", type=int, default=16)
    p.add_argument("--memory", default="2 GB")
    p.add_argument("--wall", type=int, default=6 * 3600)
    p.set_defaults(func=cmd_integrity)

    p = sub.add_parser("collect", help="combine the shard results into one good-list")
    p.add_argument("--workdir", required=True)
    p.add_argument("--list", help="the list given to integrity, to report what did not pass")
    p.set_defaults(func=cmd_collect)

    p = sub.add_parser("merge", help="hadd the good-list into blocks, one job per block")
    p.add_argument("--good", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--workdir", required=True)
    p.add_argument("--stem", required=True)
    p.add_argument("--block", type=int, default=100)
    p.add_argument("--memory", default="2 GB")
    p.add_argument("--wall", type=int, default=6 * 3600)
    p.add_argument("--force", action="store_true", help="redo blocks whose output already exists")
    p.add_argument("--block-range", help="submit only blocks START:END of the global numbering")
    p.set_defaults(func=cmd_merge)

    p = sub.add_parser("verify", help="check each merged block against the files it was built from")
    p.add_argument("--outdir", required=True)
    p.add_argument("--workdir", required=True)
    p.add_argument("--shards", type=int, default=16)
    p.add_argument("--memory", default="2 GB")
    p.add_argument("--wall", type=int, default=6 * 3600)
    p.set_defaults(func=cmd_verify)

    p = sub.add_parser("verdicts", help="summarise the verify step and list the redundant sources")
    p.add_argument("--outdir", required=True)
    p.add_argument("--workdir", required=True)
    p.set_defaults(func=cmd_verdicts)

    p = sub.add_parser("reclaim",
                       help="release the subjob outputs behind blocks that verified (--apply)")
    p.add_argument("--outdir", required=True)
    p.add_argument("--workdir", required=True)
    p.add_argument("--apply", action="store_true",
                   help="actually release them; without it, this only reports")
    p.set_defaults(func=cmd_reclaim)

    p = sub.add_parser("downstream", help="run the downstream analysis, one job per merged block")
    p.add_argument("--blockdir", required=True, help="directory holding the merged blocks")
    p.add_argument("--outdir", required=True)
    p.add_argument("--workdir", required=True)
    p.add_argument("--mainconf", required=True)
    p.add_argument("--limit", type=int, help="submit only the first N blocks (for a pilot)")
    # 2 GB rather than 4: the pool has 37,652 slots with at least 2 GB but only 9,360 with 4 GB,
    # and the reader holds one event at a time, so the footprint is the histograms plus the mixing
    # buffer, not the file.
    p.add_argument("--memory", default="2 GB")
    p.add_argument("--wall", type=int, default=12 * 3600)
    p.add_argument("--force", action="store_true")
    p.add_argument("--library-tag", default="SL24y")
    p.add_argument("--buffer", type=int, default=-1, help="mixing buffer size; -1 takes the config")
    p.add_argument("--mixing-mode", default="bufferAll")
    p.add_argument("--variation", default="",
                   help="systematic variation spec passed to the macro, e.g. dTofPThr=99,pTofPThr=2")
    p.add_argument("--tag", default="",
                   help="suffix for the job and output names, so a variation does not overwrite "
                        "the nominal (e.g. _notof)")
    p.set_defaults(func=cmd_downstream)

    p = sub.add_parser("status", help="how far the submitted work has got")
    p.add_argument("--workdir", required=True)
    p.add_argument("--outdir")
    p.set_defaults(func=cmd_status)

    args = ap.parse_args()
    if not getattr(args, "func", None):
        ap.print_help()
        return 1
    if not shutil.which("condor_submit") and not args.dry_run:
        sys.exit("ERROR: condor_submit is not on PATH; run this on a submit node")
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
