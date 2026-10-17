#!/usr/bin/env python3
"""Write a SUMS joblist whose input is a list of files instead of a catalog query.

Use it for any stage that reads files you already have: the pair stage of StCandidateFemtoMaker on candidate stores,
or a PicoDst stage on your own list (bad runs removed, a part of the production). Command, macro, mainconf, log and err
directories and provenance come from the same template and analysis info as script/generate_joblist.sh, so
job/run/submit.sh and its watch-merge tools treat the result like any generated joblist.

Differences from generate_joblist.sh:
  input        filelist:<list>, fileListSyntax "paths" for local files, "xrootd" for root:// URLs
  per subjob   at most --files-per-subjob input files (SUMS spreads the files evenly over the subjobs)
  outputs      <stem>_$JOBID.root to --output-dir; with --tree-output-dir also <stem>_$JOBID.tree.root there
  sandbox      zipped at submission (installer="ZIP", a new package name each time this script runs), so later edits of
               the working tree do not reach subjobs that are still idle
  file name    job/joblist/joblist_<anaName>_<label>.xml (job/run/submit.sh takes <anaName>_<label> as its name)

--input is a list file (one path or root:// URL per line; '#' starts a comment) or a directory. A directory gives every
*.root file in it with size > 0, written next to the outputs as <output-dir>/<anaName>_<label>_input_files.list.

Examples for StCandidateFemtoMaker (repository root; <v> = rootfile/<anaName>/<major>):
  # 1st stage: candidate stores from a PicoDst list
  ./script/generate_joblist_from_list.py --mainconf config/mainconf/main_<anaName>_candidateStore.yaml \\
      --input <PicoDst list> --files-per-subjob 3 --label candidateStore --output-dir <v>/<version>/store
  # pair stage on the stores: candidate trees for BDT training
  ./script/generate_joblist_from_list.py --mainconf config/mainconf/main_<anaName>_bdtTraining.yaml \\
      --input <v>/<store version>/store --files-per-subjob 100 --label bdtTrainingTrees \\
      --output-dir <v>/<version>/hist --tree-output-dir <v>/<version>/tree
  # pair stage on the stores: production with BDT scores
  ./script/generate_joblist_from_list.py --mainconf config/mainconf/main_<anaName>_bdt.yaml \\
      --input <v>/<store version>/store --files-per-subjob 100 --label bdtAnalysis --output-dir <v>/<version>/hist
Then: cd job/run && ./submit.sh ../joblist/joblist_<anaName>_<label>.xml
Regenerate right before submitting, so that the recorded git revision describes the submitted tree.
"""
import argparse
import datetime
import glob
import os
import re
import subprocess
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
import analysis_info_helper as helper  # noqa: E402

PROJECT_ROOT = helper.get_project_root(SCRIPT_DIR)


def fail(message):
    sys.exit("ERROR: " + message)


def project_path(path):
    """Absolute path; relative paths are relative to the repository root."""
    path = os.path.expandvars(os.path.expanduser(path))
    return os.path.abspath(path if os.path.isabs(path) else os.path.join(PROJECT_ROOT, path))


def git(args):
    try:
        return subprocess.check_output(['git', '-C', PROJECT_ROOT] + args, stderr=subprocess.DEVNULL, text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        return ''


def read_input(path, list_out):
    """Return (list file to give SUMS, entries). A directory is turned into a list of its non-empty *.root files."""
    if os.path.isdir(path):
        files = sorted(glob.glob(os.path.join(path, '*.root')))
        entries = [f for f in files if os.path.getsize(f) > 0]
        if len(files) > len(entries):
            print("skipped {} empty files (unfinished subjobs?) in {}".format(len(files) - len(entries), path))
        if not entries:
            fail("no non-empty *.root files in " + path)
        os.makedirs(os.path.dirname(list_out), exist_ok=True)
        with open(list_out, 'w') as f:
            f.write('\n'.join(entries) + '\n')
        return list_out, entries
    if not os.path.isfile(path):
        fail("input is neither a list file nor a directory: " + path)
    with open(path) as f:
        entries = [line.split('#')[0].split()[0] for line in f if line.split('#')[0].strip()]
    if not entries:
        fail("the list is empty: " + path)
    return path, entries


def file_list_syntax(entries):
    remote = [e for e in entries if e.startswith('root://')]
    if remote and len(remote) != len(entries):
        fail("the list mixes root:// URLs and local paths")
    if remote:
        return 'xrootd'
    if any(not os.path.isabs(e) for e in entries):
        fail("local paths in the list must be absolute (the farm node does not run in this directory)")
    print("checking that the {} listed files exist and are not empty ...".format(len(entries)))
    bad = [e for e in entries if not os.path.isfile(e) or os.path.getsize(e) == 0]
    if bad:
        fail("{} of {} listed files are missing or empty, e.g. {}".format(len(bad), len(entries), bad[0]))
    return 'paths'


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0], epilog=__doc__.split('\n', 2)[2],
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--mainconf', required=True, help='config/mainconf/<...>.yaml run by every subjob')
    parser.add_argument('--input', required=True, help='list file or directory of input files')
    parser.add_argument('--files-per-subjob', required=True, type=int, help='at most this many input files per subjob (SUMS balances them)')
    parser.add_argument('--label', required=True, help='what this job is, e.g. bdtAnalysis (letters and digits)')
    parser.add_argument('--output-dir', required=True, help='directory for <stem>_$JOBID.root')
    parser.add_argument('--tree-output-dir', help='directory for <stem>_$JOBID.tree.root (jobs writing candidate trees)')
    args = parser.parse_args()

    if not re.match(r'^[A-Za-z0-9]+$', args.label):
        fail("--label must be letters and digits only: " + args.label)
    if args.files_per_subjob < 1:
        fail("--files-per-subjob must be at least 1")

    mainconf_path = project_path(args.mainconf)
    if not os.path.isfile(mainconf_path):
        fail("mainconf not found: " + mainconf_path)
    mainconf = helper.get_mainconf_for_command(helper.normalize_relpath(os.path.relpath(mainconf_path, PROJECT_ROOT)))
    if not mainconf.startswith('config/mainconf/'):
        fail("the mainconf must be under config/mainconf/ (submit.sh reads it back from the joblist): " + mainconf)
    analysis_path = helper.get_analysis_path_from_mainconf(mainconf_path, os.path.join(PROJECT_ROOT, 'config'))
    if not os.path.isfile(analysis_path):
        fail("analysis info not found: " + analysis_path)
    info = helper.load_analysis_info(analysis_path)
    analysis = info.get('analysis') or {}
    star_tag = info.get('starTag') or info.get('startag') or {}
    ana_name = analysis.get('anaName') or analysis.get('name')
    for key, value in (('analysis.anaName', ana_name), ('analysis.baseRunMacro', analysis.get('baseRunMacro')),
                       ('analysis.baseAnaMacro', analysis.get('baseAnaMacro')), ('starTag.libraryTag', star_tag.get('libraryTag'))):
        if not value:
            fail("{} is not set in {}".format(key, analysis_path))

    joblist_dir = os.path.join(PROJECT_ROOT, 'job', 'joblist')
    output_dir = project_path(args.output_dir)
    tree_dir = project_path(args.tree_output_dir) if args.tree_output_dir else None
    list_file, entries = read_input(project_path(args.input),
                                    os.path.join(output_dir, '{}_{}_input_files.list'.format(ana_name, args.label)))
    syntax = file_list_syntax(entries)
    for d in filter(None, (output_dir, tree_dir)):
        if glob.glob(os.path.join(d, '*.root')):
            print("note: {} already holds ROOT files; the new outputs are told apart by their JOBID".format(d))

    work_dir = helper.resolve_work_dir(analysis, PROJECT_ROOT)
    stem = analysis.get('outputFileStem') or ana_name
    scratch_subdir = analysis.get('scratchSubdir') or ana_name
    max_events = analysis.get('maxEvents', -1)
    if max_events is None or str(max_events).strip() == '':
        max_events = -1
    rev = git(['rev-parse', 'HEAD']) or 'unknown'
    dirty = str(len([l for l in git(['status', '--porcelain']).splitlines() if not l.startswith('??')]))

    # Same placeholders as analysis_info_helper.py --generate-joblist; the input placeholders are replaced below.
    with open(os.path.join(joblist_dir, 'job_template_from_conf.xml')) as f:
        content = f.read()
    for placeholder, value in (
            ('__GIT_REV__', rev), ('__GIT_DIRTY__', dirty),
            ('__JOB_NAME__', '{}_{}'.format(analysis.get('jobName') or ana_name, args.label)),
            ('__RUN_MACRO__', analysis['baseRunMacro'] + '.C'), ('__STARVER__', star_tag['libraryTag']),
            ('__SCRATCH_SUBDIR__', scratch_subdir), ('__OUTPUT_FILE_STEM__', stem), ('__MAINCONF__', mainconf),
            ('__MAX_EVENTS__', str(max_events)), ('__WORK_DIR__', work_dir),
            ('__LOG_DIR__', helper.resolve_batch_dir(analysis, 'logDir', os.path.join(work_dir, 'log'), PROJECT_ROOT)),
            ('__ERR_DIR__', helper.resolve_batch_dir(analysis, 'errDir', os.path.join(work_dir, 'err'), PROJECT_ROOT)),
            ('__ANA_SO_PREFIX__', analysis['baseAnaMacro']), ('__FILES_PER_JOB__', str(args.files_per_subjob))):
        content = content.replace(placeholder, value)

    package = '{}_{}_{}'.format(ana_name, args.label, datetime.datetime.now().strftime('%Y%m%d%H%M%S'))
    output_line = '<output fromScratch="{0}/{1}_$JOBID{2}" toURL="file:{3}/" />'
    outputs = [output_line.format(scratch_subdir, stem, '.root', output_dir)]
    if tree_dir:
        outputs.append(output_line.format(scratch_subdir, stem, '.tree.root', tree_dir))
    comment = ('<!--\n  Generated by script/generate_joblist_from_list.py from job/joblist/job_template_from_conf.xml:\n'
               '  {} input files from {}, {} per subjob, mainconf {}.\n'
               '  Sandbox zipped at submission (package {}).\n-->').format(
                   len(entries), os.path.relpath(list_file, PROJECT_ROOT), args.files_per_subjob, mainconf, package)
    edits = (
        (r'<!--.*?-->', comment),
        (r'fileListSyntax="[a-z]+"', 'fileListSyntax="{}"'.format(syntax)),
        (r'<input\s[^>]*/>', '<input URL="filelist:{}" />'.format(list_file)),
        (r'<SandBox>', '<SandBox installer="ZIP">'),
        (r'<Package>', '<Package name="{}">'.format(package)),
        (r'<output\s[^>]*/>', '\n  '.join(outputs)),
    )
    for pattern, replacement in edits:
        content, n = re.subn(pattern, lambda m: replacement, content, count=1, flags=re.S)
        if n != 1:
            fail("the template job/joblist/job_template_from_conf.xml no longer matches " + pattern)
    left = sorted(set(re.findall(r'__[A-Z_]+__', content)))
    if left:
        fail("the template has placeholders this script does not fill: " + ", ".join(left))

    out = os.path.join(joblist_dir, 'joblist_{}_{}.xml'.format(ana_name, args.label))
    with open(out, 'w') as f:
        f.write(content)
    print("wrote {}: {} input files ({}), {} per subjob, about {} subjobs".format(
        os.path.relpath(out, PROJECT_ROOT), len(entries), syntax, args.files_per_subjob,
        -(-len(entries) // args.files_per_subjob)))
    print("  mainconf {}\n  outputs  {}{}".format(mainconf, output_dir, '\n           ' + tree_dir if tree_dir else ''))
    print("submit: cd job/run && ./submit.sh ../joblist/{}".format(os.path.basename(out)))


if __name__ == '__main__':
    main()
