#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Read mainconf and its analysis info YAML; output libraryTag for setup.sh,
generate joblist XML, or extract the embedded mainconf from a joblist.
Python 2.7 compatible.
Usage:
  python script/analysis_info_helper.py --library-tag [--mainconf PATH]
  python script/analysis_info_helper.py --generate-joblist MAINCONF_PATH [--project-root PATH]
  python script/analysis_info_helper.py --mainconf-from-joblist JOBLIST_XML [--project-root PATH]
Example:
  python script/analysis_info_helper.py --generate-joblist config/mainconf/main_auau19_anaLambda.yaml
"""
from __future__ import print_function
import os
import re
import sys
import argparse

try:
    import yaml
except ImportError:
    yaml = None


def get_project_root(script_dir):
    """Script is in script/; project root is parent of script/."""
    return os.path.dirname(script_dir)


def load_yaml(path):
    if yaml is None:
        print("ERROR: PyYAML is required. Install with: pip install pyyaml", file=sys.stderr)
        sys.exit(1)
    with open(path, 'r') as f:
        return yaml.safe_load(f)


def get_analysis_path_from_mainconf_no_yaml(mainconf_path, config_base):
    """When PyYAML is missing: grep for 'analysis:' in mainconf (value relative to config/)."""
    with open(mainconf_path, 'r') as f:
        for line in f:
            m = re.search(r'analysis\s*:\s*(\S+)', line)
            if m:
                return os.path.join(config_base, m.group(1).strip())
    return None


def get_library_tag_from_file_no_yaml(analysis_path):
    """When PyYAML is missing: grep for 'libraryTag:' in analysis file (under starTag: block)."""
    with open(analysis_path, 'r') as f:
        for line in f:
            m = re.search(r'libraryTag\s*:\s*(\S+)', line)
            if m:
                return m.group(1).strip()
    return None


def get_ana_name_from_file_no_yaml(analysis_path):
    """When PyYAML is missing: grep for 'anaName:' and extract quoted value."""
    with open(analysis_path, 'r') as f:
        for line in f:
            if 'anaName' not in line:
                continue
            m = re.search(r'anaName\s*:\s*(?:&\w+\s+)?["\']([^"\']+)["\']', line)
            if m:
                return m.group(1).strip()
            m = re.search(r'anaName\s*:\s*(\S+)', line)
            if m:
                return m.group(1).strip()
    return None


def get_all_pico_dst_list_no_yaml(analysis_path):
    """When PyYAML is missing: grep for 'allPicoDstList:' in analysis file."""
    with open(analysis_path, 'r') as f:
        for line in f:
            m = re.search(r'allPicoDstList\s*:\s*["\']?([^"\'\s]+)["\']?', line)
            if m:
                return m.group(1).strip()
    return None


def get_analysis_path_from_mainconf(mainconf_path, config_base):
    """mainconf has 'analysis: path'; path is relative to config/."""
    data = load_yaml(mainconf_path)
    rel = data.get('analysis') or data.get('analysis info')
    if not rel:
        print("ERROR: mainconf has no 'analysis' key: {}".format(mainconf_path), file=sys.stderr)
        sys.exit(1)
    rel = rel.strip()
    return os.path.join(config_base, rel)


def load_analysis_info(analysis_path):
    return load_yaml(analysis_path)


def normalize_relpath(path):
    return path.replace(os.sep, '/')


def get_mainconf_relpath(mainconf_arg, mainconf_path, project_root):
    """Return project-relative mainconf path using forward slashes."""
    if os.path.isabs(mainconf_arg):
        rel = os.path.relpath(mainconf_path, project_root)
    else:
        rel = mainconf_arg
    return normalize_relpath(rel)


def get_mainconf_for_command(mainconf_rel):
    """Return config/... mainconf path for joblist command."""
    if mainconf_rel.startswith('config/'):
        return mainconf_rel
    return 'config/' + mainconf_rel


def extract_mainconf_from_joblist(joblist_path):
    """Extract the embedded config/mainconf/...yaml path from a joblist XML."""
    with open(joblist_path, 'r') as f:
        content = f.read()
    matches = re.findall(r'(config/mainconf/[A-Za-z0-9_./-]+\.ya?ml)', content)
    unique = []
    for match in matches:
        if match not in unique:
            unique.append(match)
    if not unique:
        print("ERROR: embedded mainconf not found in joblist: {}".format(joblist_path), file=sys.stderr)
        sys.exit(1)
    if len(unique) > 1:
        print("ERROR: multiple embedded mainconf paths found in joblist: {}".format(", ".join(unique)), file=sys.stderr)
        sys.exit(1)
    return unique[0]


def build_catalog_url(star_tag):
    """Build SUMS catalog URL from starTag."""
    base = "catalog:star.bnl.gov"
    parts = []
    if star_tag.get('triggerSets'):
        parts.append("trgsetupname={}".format(star_tag['triggerSets']))
    if star_tag.get('productionTag'):
        parts.append("production={}".format(star_tag['productionTag']))
    if star_tag.get('filetype'):
        parts.append("filetype={}".format(star_tag['filetype']))
    if star_tag.get('filenameFilter'):
        parts.append("filename~{}".format(star_tag['filenameFilter']))
    if star_tag.get('filenameNotFilter'):
        parts.append("filename!~{}".format(star_tag['filenameNotFilter']))
    if star_tag.get('storageExclude'):
        parts.append("storage!={}".format(star_tag['storageExclude']))
    if star_tag.get('storage'):
        parts.append("storage={}".format(star_tag['storage']))
    return base + "?" + ",".join(parts)


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = get_project_root(script_dir)
    config_base = os.path.join(project_root, 'config')

    parser = argparse.ArgumentParser(description='Read analysis info from mainconf; output libraryTag, generate joblist, or inspect joblist mainconf.')
    parser.add_argument('mainconf', nargs='?', default=None,
                        help='Main config path (e.g. config/mainconf/main_auau19_anaLambda.yaml). For --generate-joblist, pass as argument; for --library-tag, optional (else uses --mainconf).')
    parser.add_argument('--mainconf', dest='mainconf_opt', default=None,
                        help='Main config path (used by setup.sh). Ignored if mainconf is passed as positional.')
    parser.add_argument('--project-root', default=project_root, help='Project root directory')
    parser.add_argument('--library-tag', action='store_true', help='Print starTag.libraryTag to stdout')
    parser.add_argument('--ana-name', action='store_true', help='Print analysis anaName from analysis_info to stdout')
    parser.add_argument('--pico-dst-list', action='store_true', help='Print default picoDst list path (config/...) for run script')
    parser.add_argument('--output-rootfile', action='store_true', help='Print default output root path: rootfile/anaName/anaName_temp.root')
    parser.add_argument('--generate-joblist', action='store_true', help='Generate joblist XML from template')
    parser.add_argument('--mainconf-from-joblist', dest='joblist_path', default=None,
                        help='Extract embedded config/mainconf/...yaml path from a joblist XML')
    args = parser.parse_args()

    project_root = os.path.abspath(args.project_root)
    config_base = os.path.join(project_root, 'config')

    if args.joblist_path:
        joblist_path = args.joblist_path
        if not os.path.isabs(joblist_path):
            joblist_path = os.path.join(project_root, joblist_path)
        if not os.path.isfile(joblist_path):
            print("ERROR: joblist not found: {}".format(joblist_path), file=sys.stderr)
            sys.exit(1)
        print(extract_mainconf_from_joblist(joblist_path))
        return

    mainconf_arg = args.mainconf
    if mainconf_arg is None:
        mainconf_arg = args.mainconf_opt
    if mainconf_arg is None and not args.generate_joblist:
        mainconf_arg = 'config/mainconf/main_auau19_anaLambda.yaml'
    if args.generate_joblist and (mainconf_arg is None or not mainconf_arg.strip()):
        print("ERROR: --generate-joblist requires mainconf. Pass it as the first argument or use --mainconf.", file=sys.stderr)
        sys.exit(1)
    mainconf_path = mainconf_arg
    if not os.path.isabs(mainconf_path):
        mainconf_path = os.path.join(project_root, mainconf_path)
    if not os.path.isfile(mainconf_path):
        print("ERROR: mainconf not found: {}".format(mainconf_path), file=sys.stderr)
        sys.exit(1)

    if yaml is not None:
        analysis_path = get_analysis_path_from_mainconf(mainconf_path, config_base)
        if not os.path.isfile(analysis_path):
            print("ERROR: analysis info not found: {}".format(analysis_path), file=sys.stderr)
            sys.exit(1)
        info = load_analysis_info(analysis_path)
        star_tag = info.get('starTag') or info.get('startag') or {}
        analysis = info.get('analysis') or {}
    else:
        analysis_path = get_analysis_path_from_mainconf_no_yaml(mainconf_path, config_base)
        if not analysis_path or not os.path.isfile(analysis_path):
            print("ERROR: analysis info not found (mainconf: {})".format(mainconf_path), file=sys.stderr)
            sys.exit(1)
        star_tag = {}
        analysis = {}

    if args.library_tag:
        if yaml is not None:
            tag = star_tag.get('libraryTag', 'SL24y')
        else:
            tag = get_library_tag_from_file_no_yaml(analysis_path)
            if tag is None:
                tag = 'SL24y'
        print(tag)
        return

    if args.ana_name:
        if yaml is not None:
            name = analysis.get('anaName') or analysis.get('name') or 'auau19_anaPhi'
        else:
            name = get_ana_name_from_file_no_yaml(analysis_path)
            if name is None:
                name = 'auau19_anaPhi'
        print(name)
        return

    if args.pico_dst_list:
        if yaml is not None:
            dataset = info.get('dataset') or {}
            list_rel = dataset.get('allPicoDstList', 'picoDstList/auau19GeV.list')
        else:
            list_rel = get_all_pico_dst_list_no_yaml(analysis_path)
            if list_rel is None:
                list_rel = 'picoDstList/auau19GeV.list'
        print(os.path.join('config', list_rel))
        return

    if args.output_rootfile:
        if yaml is not None:
            name = analysis.get('anaName') or analysis.get('name') or 'auau19_anaPhi'
        else:
            name = get_ana_name_from_file_no_yaml(analysis_path)
            if name is None:
                name = 'auau19_anaPhi'
        print(os.path.join('rootfile', name, name + '_temp.root'))
        return

    if args.generate_joblist:
        if yaml is None:
            print("ERROR: --generate-joblist requires PyYAML. Install with: pip install pyyaml", file=sys.stderr)
            sys.exit(1)
        mainconf_rel = get_mainconf_relpath(mainconf_arg, mainconf_path, project_root)
        mainconf_for_command = get_mainconf_for_command(mainconf_rel)

        base_run = analysis.get('baseRunMacro', 'run_anaLambda')
        base_ana = analysis.get('baseAnaMacro', 'anaLambda')
        run_macro = base_run + '.C'
        ana_so_prefix = base_ana

        ana_name = analysis.get('anaName') or analysis.get('name') or 'auau19_anaLambda_temp'
        job_name = analysis.get('jobName') or ana_name
        scratch_subdir = analysis.get('scratchSubdir') or ana_name
        output_stem = analysis.get('outputFileStem') or ana_name
        n_files = analysis.get('nFiles', 1000)
        max_events = analysis.get('maxEvents', -1)
        if max_events is None or str(max_events).strip() == '':
            max_events = -1
        work_dir = analysis.get('workDir', '/star/u/$USER/Path/To/star-analyzer')
        starver = star_tag.get('libraryTag', 'SL24y')
        catalog_url = build_catalog_url(star_tag)

        template_path = os.path.join(project_root, 'job', 'joblist', 'job_template_from_conf.xml')
        if not os.path.isfile(template_path):
            print("ERROR: template not found: {}".format(template_path), file=sys.stderr)
            sys.exit(1)
        with open(template_path, 'r') as f:
            content = f.read()

        replacements = [
            ('__JOB_NAME__', job_name),
            ('__RUN_MACRO__', run_macro),
            ('__STARVER__', starver),
            ('__SCRATCH_SUBDIR__', scratch_subdir),
            ('__OUTPUT_FILE_STEM__', output_stem),
            ('__MAINCONF__', mainconf_for_command),
            ('__MAX_EVENTS__', str(max_events)),
            ('__WORK_DIR__', work_dir),
            ('__CATALOG_URL__', catalog_url),
            ('__N_FILES__', str(n_files)),
            ('__ANA_SO_PREFIX__', ana_so_prefix),
        ]
        for placeholder, value in replacements:
            content = content.replace(placeholder, value)

        out_basename = 'joblist_' + ana_name + '.xml'
        out_dir = os.path.join(project_root, 'job', 'joblist')
        out_path = os.path.join(out_dir, out_basename)
        with open(out_path, 'w') as f:
            f.write(content)
        print("Wrote {}".format(out_path))
        return

    print("ERROR: specify --library-tag, --ana-name, --pico-dst-list, --output-rootfile, or --generate-joblist", file=sys.stderr)
    sys.exit(1)


if __name__ == '__main__':
    main()
