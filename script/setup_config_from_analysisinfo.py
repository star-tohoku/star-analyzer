#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Create analysis-specific config files from templates using anaName from analysis_info.
Usage:
  python script/setup_config_from_analysisinfo.py ANALYSIS_INFO [--config-dir PATH]
Example:
  python script/setup_config_from_analysisinfo.py config/analysis/analysis_info_pp500_anaPhi.yaml

Reads analysis_info YAML to get anaName, then:
Copies:
  cuts/event/event.yaml       -> cuts/event/event_{anaName}.yaml
  cuts/mixing/mixing.yaml     -> cuts/mixing/mixing_{anaName}.yaml
  cuts/pid/pid.yaml           -> cuts/pid/pid_{anaName}.yaml
  cuts/track/track.yaml       -> cuts/track/track_{anaName}.yaml
  cuts/v0reco/v0.yaml         -> cuts/v0reco/v0_{anaName}.yaml
  maker/maker_anaPhi.yaml     -> maker/maker_{anaName}.yaml
Creates:
  mainconf/main_{anaName}.yaml  (from mainconf/mainconf.yaml; __ANANAME__ placeholder is replaced)
If a destination already exists, it is overwritten.
"""
from __future__ import print_function
import os
import re
import shutil
import sys
import argparse

try:
    import yaml
except ImportError:
    yaml = None


def get_project_root(script_dir):
    """Project root is parent of script/."""
    return os.path.dirname(script_dir)


def load_yaml(path):
    if yaml is None:
        print("ERROR: PyYAML is required. Install with: pip install pyyaml", file=sys.stderr)
        sys.exit(1)
    with open(path, 'r') as f:
        return yaml.safe_load(f)


def get_analysis_field_no_yaml(analysis_path, key):
    """When PyYAML is missing: grep for key and extract value."""
    with open(analysis_path, 'r') as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line or line.startswith('#'):
                continue
            # Ignore trailing inline comments before matching key/value.
            line = line.split('#', 1)[0].strip()
            if key not in line:
                continue
            m = re.search(r'{}\s*:\s*(?:&\w+\s+)?["\']([^"\']+)["\']'.format(key), line)
            if m:
                return m.group(1).strip()
            m = re.search(r'{}\s*:\s*(\S+)'.format(key), line)
            if m:
                return m.group(1).strip()
    return None


# (src relative to config/, dst subdir under config/, base filename for dst)
CUTS_COPIES = [
    ('cuts/event/event.yaml', 'cuts/event', 'event'),
    ('cuts/mixing/mixing.yaml', 'cuts/mixing', 'mixing'),
    ('cuts/pid/pid.yaml', 'cuts/pid', 'pid'),
    ('cuts/track/track.yaml', 'cuts/track', 'track'),
    ('cuts/v0reco/v0.yaml', 'cuts/v0reco', 'v0'),
    ('cuts/centrality/centrality.yaml', 'cuts/centrality', 'centrality'),
]
MAKER_SRC = 'maker/maker_anaPhi.yaml'
MAKER_DIR = 'maker'
MAKER_BASE = 'maker'
MAINCONF_TEMPLATE_REL = 'mainconf/mainconf.yaml'
MAINCONF_TEMPLATE_ANANAME = '__ANANAME__'


def write_mainconf(config_base, ana_name, analysis_rel):
    """Create mainconf/main_{anaName}.yaml from template, with analysis: -> analysis_rel."""
    template_path = os.path.join(config_base, MAINCONF_TEMPLATE_REL)
    mainconf_dir = os.path.join(config_base, 'mainconf')
    mainconf_dst = os.path.join(mainconf_dir, 'main_{}.yaml'.format(ana_name))
    if not os.path.isfile(template_path):
        print("WARNING: mainconf template not found, skip: {}".format(MAINCONF_TEMPLATE_REL))
        return
    with open(template_path, 'r') as f:
        content = f.read()
    content = content.replace(MAINCONF_TEMPLATE_ANANAME, ana_name)
    
    # Determine the maker key based on analysis name (e.g. lambda:, nuclearid:, phi:)
    maker_key = 'lambda'
    if 'nuclearid' in ana_name.lower():
        maker_key = 'nuclearid'
    elif 'phi' in ana_name.lower():
        maker_key = 'phi'
    elif 'lambda' in ana_name.lower():
        maker_key = 'lambda'

    content = re.sub(
        r'^lambda\s*:\s*(maker/maker_{}\.yaml)'.format(re.escape(ana_name)),
        '{}:        \\1'.format(maker_key),
        content,
        flags=re.MULTILINE
    )

    # Point to analysis-specific cut configs we create
    for _subdir, _base in [('event', 'event'), ('track', 'track'), ('pid', 'pid'),
                           ('v0reco', 'v0'), ('mixing', 'mixing'), ('centrality', 'centrality')]:
        content = content.replace(
            'cuts/{}/{}.yaml'.format(_subdir, _base),
            'cuts/{}/{}_{}.yaml'.format(_subdir, _base, ana_name)
        )
    content = re.sub(
        r'^analysis\s*:\s*\S+',
        'analysis:      {}'.format(analysis_rel),
        content,
        flags=re.MULTILINE
    )
    if not os.path.isdir(mainconf_dir):
        try:
            os.makedirs(mainconf_dir)
        except OSError:
            if not os.path.isdir(mainconf_dir):
                raise
    with open(mainconf_dst, 'w') as f:
        f.write(content)
    print("Created: {}".format(os.path.relpath(mainconf_dst, config_base)))


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = get_project_root(script_dir)
    default_config = os.path.join(project_root, 'config')

    parser = argparse.ArgumentParser(
        description='Create analysis-specific config files from analysis_info (using anaName).'
    )
    parser.add_argument(
        'analysis_info',
        help='Analysis info YAML path (e.g. config/analysis/analysis_info_pp500_anaPhi.yaml)'
    )
    parser.add_argument(
        '--config-dir',
        default=default_config,
        help='Config directory (default: project_root/config)'
    )
    args = parser.parse_args()

    config_base = os.path.abspath(args.config_dir)
    analysis_path = args.analysis_info
    if not os.path.isabs(analysis_path):
        analysis_path = os.path.join(project_root, analysis_path)
    if not os.path.isfile(analysis_path):
        print("ERROR: analysis_info not found: {}".format(analysis_path), file=sys.stderr)
        sys.exit(1)
    analysis_rel = os.path.relpath(analysis_path, config_base).replace(os.sep, '/')

    if yaml is not None:
        info = load_yaml(analysis_path)
        analysis = info.get('analysis') or {}
        ana_name = analysis.get('anaName') or analysis.get('name')
        base_ana_macro = analysis.get('baseAnaMacro')
        base_run_macro = analysis.get('baseRunMacro')
    else:
        ana_name = get_analysis_field_no_yaml(analysis_path, 'anaName')
        base_ana_macro = get_analysis_field_no_yaml(analysis_path, 'baseAnaMacro')
        base_run_macro = get_analysis_field_no_yaml(analysis_path, 'baseRunMacro')
    if not ana_name:
        print("ERROR: anaName not found in analysis info: {}".format(analysis_path), file=sys.stderr)
        sys.exit(1)

    print("anaName: {}".format(ana_name))
    print("config base: {}".format(config_base))

    def copy_one(src_rel, subdir, base_name):
        src = os.path.join(config_base, src_rel)
        dst_name = "{}_{}.yaml".format(base_name, ana_name)
        dst = os.path.join(config_base, subdir, dst_name)
        if not os.path.isfile(src):
            print("WARNING: template not found, skip: {}".format(src_rel))
            return
        d = os.path.dirname(dst)
        if not os.path.isdir(d):
            try:
                os.makedirs(d)
            except OSError:
                if not os.path.isdir(d):
                    raise
        shutil.copy2(src, dst)
        print("Created: {}".format(os.path.relpath(dst, config_base)))

    for src_rel, subdir, base_name in CUTS_COPIES:
        copy_one(src_rel, subdir, base_name)

    # maker: template based on analysis name -> maker/maker_{anaName}.yaml
    maker_template_rel = 'maker/maker_anaPhi.yaml'
    if 'nuclearid' in ana_name.lower():
        maker_template_rel = 'maker/maker_anaNuclearId.yaml'
    elif 'lambda' in ana_name.lower():
        maker_template_rel = 'maker/maker_anaLambda.yaml'
    elif 'phi' in ana_name.lower():
        maker_template_rel = 'maker/maker_anaPhi.yaml'

    maker_src = os.path.join(config_base, maker_template_rel)
    maker_dst = os.path.join(config_base, MAKER_DIR, "{}_{}.yaml".format(MAKER_BASE, ana_name))
    if not os.path.isfile(maker_src):
        print("WARNING: template not found, skip: {}".format(maker_template_rel))
    else:
        shutil.copy2(maker_src, maker_dst)
        print("Created: {}".format(os.path.relpath(maker_dst, config_base)))

    write_mainconf(config_base, ana_name, analysis_rel)

    # generate macros if base macros are defined
    if base_ana_macro and base_run_macro:
        maker_name = base_ana_macro[3:] if base_ana_macro.startswith('ana') else base_ana_macro
        maker_class = "St{}Maker".format(maker_name)
        maker_var = "{}Maker".format(maker_name.lower())
        
        # ana macro
        ana_src = os.path.join(project_root, 'analysis', 'anaPhi.C')
        ana_dst = os.path.join(project_root, 'analysis', '{}.C'.format(base_ana_macro))
        if os.path.isfile(ana_src):
            with open(ana_src, 'r') as f:
                content = f.read()
            content = content.replace('anaPhi', base_ana_macro)
            content = content.replace('StPhiMaker', maker_class)
            content = content.replace('phiMaker', maker_var)
            content = content.replace('"phi"', '"{}"'.format(maker_name.lower()))
            content = content.replace('auau19_anaPhi', ana_name)
            content = content.replace('run_anaPhi', base_run_macro)
            with open(ana_dst, 'w') as f:
                f.write(content)
            print("Created: {}".format(os.path.relpath(ana_dst, project_root)))
            
        # run macro
        run_src = os.path.join(project_root, 'analysis', 'run_anaPhi.C')
        run_dst = os.path.join(project_root, 'analysis', '{}.C'.format(base_run_macro))
        if os.path.isfile(run_src):
            with open(run_src, 'r') as f:
                content = f.read()
            content = content.replace('run_anaPhi', base_run_macro)
            content = content.replace('anaPhi', base_ana_macro)
            content = content.replace('StPhiMaker', maker_class)
            content = content.replace('auau19_anaPhi', ana_name)
            with open(run_dst, 'w') as f:
                f.write(content)
            print("Created: {}".format(os.path.relpath(run_dst, project_root)))


if __name__ == '__main__':
    main()
