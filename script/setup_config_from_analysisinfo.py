#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Create analysis-specific config files from templates using anaName from analysis_info.
Usage:
  python script/setup_config_from_mainconf.py ANALYSIS_INFO [--config-dir PATH]
Example:
  python script/setup_config_from_mainconf.py config/analysis/analysis_info_pp500_anaPhi.yaml

Reads analysis_info YAML to get anaName, then:
Copies:
  cuts/event/event.yaml       -> cuts/event/event_{anaName}.yaml
  cuts/mixing/mixing.yaml     -> cuts/mixing/mixing_{anaName}.yaml
  cuts/pid/pid.yaml           -> cuts/pid/pid_{anaName}.yaml
  cuts/track/track.yaml       -> cuts/track/track_{anaName}.yaml
  cuts/v0reco/v0.yaml         -> cuts/v0reco/v0_{anaName}.yaml
  maker/maker_anaPhi.yaml     -> maker/maker_{anaName}.yaml
Creates:
  mainconf/main_{anaName}.yaml  (paths point to the new config files and given analysis_info)
If a destination already exists, prints a warning and skips that file.
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


def get_ana_name_no_yaml(analysis_path):
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


# (src relative to config/, dst subdir under config/, base filename for dst)
CUTS_COPIES = [
    ('cuts/event/event.yaml', 'cuts/event', 'event'),
    ('cuts/mixing/mixing.yaml', 'cuts/mixing', 'mixing'),
    ('cuts/pid/pid.yaml', 'cuts/pid', 'pid'),
    ('cuts/track/track.yaml', 'cuts/track', 'track'),
    ('cuts/v0reco/v0.yaml', 'cuts/v0reco', 'v0'),
]
MAKER_SRC = 'maker/maker_anaPhi.yaml'
MAKER_DIR = 'maker'
MAKER_BASE = 'maker'
MAINCONF_TEMPLATE_REL = 'mainconf/main_auau19_anaPhi.yaml'
MAINCONF_TEMPLATE_ANANAME = 'auau19_anaPhi'


def write_mainconf(config_base, ana_name, analysis_rel):
    """Create mainconf/main_{anaName}.yaml from template, with analysis: -> analysis_rel."""
    template_path = os.path.join(config_base, MAINCONF_TEMPLATE_REL)
    mainconf_dir = os.path.join(config_base, 'mainconf')
    mainconf_dst = os.path.join(mainconf_dir, 'main_{}.yaml'.format(ana_name))
    if os.path.isfile(mainconf_dst):
        print("WARNING: already exists, skip: {}".format(os.path.relpath(mainconf_dst, config_base)))
        return
    if not os.path.isfile(template_path):
        print("WARNING: mainconf template not found, skip: {}".format(MAINCONF_TEMPLATE_REL))
        return
    with open(template_path, 'r') as f:
        content = f.read()
    content = content.replace(MAINCONF_TEMPLATE_ANANAME, ana_name)
    # Point to analysis-specific cut configs we create
    for _subdir, _base in [('event', 'event'), ('track', 'track'), ('pid', 'pid'),
                           ('v0reco', 'v0'), ('mixing', 'mixing')]:
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
    else:
        ana_name = get_ana_name_no_yaml(analysis_path)
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
        if os.path.isfile(dst):
            print("WARNING: already exists, skip: {}".format(os.path.relpath(dst, config_base)))
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

    # maker: maker/maker_anaPhi.yaml -> maker/maker_{anaName}.yaml
    maker_src = os.path.join(config_base, MAKER_SRC)
    maker_dst = os.path.join(config_base, MAKER_DIR, "{}_{}.yaml".format(MAKER_BASE, ana_name))
    if not os.path.isfile(maker_src):
        print("WARNING: template not found, skip: {}".format(MAKER_SRC))
    elif os.path.isfile(maker_dst):
        print("WARNING: already exists, skip: {}".format(os.path.relpath(maker_dst, config_base)))
    else:
        shutil.copy2(maker_src, maker_dst)
        print("Created: {}".format(os.path.relpath(maker_dst, config_base)))

    write_mainconf(config_base, ana_name, analysis_rel)


if __name__ == '__main__':
    main()
