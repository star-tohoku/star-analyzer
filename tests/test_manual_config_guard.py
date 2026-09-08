#!/usr/bin/env python3
"""Isolated regression tests: never generate configs/macros in the real repo.

Run: python3 -B tests/test_manual_config_guard.py
All writes and generator subprocesses are confined to TemporaryDirectory trees.
"""
import os
from pathlib import Path
import runpy
import shutil
import subprocess
import sys
import tempfile
import unittest

sys.dont_write_bytecode = True
GENERATOR = Path(__file__).resolve().parents[1] / "script" / "setup_config_from_analysisinfo.py"
MARKER = "# config-generator: manual"


class Fixture:
    def __init__(self, root, main_text):
        self.root = Path(root)
        self.config = self.root / "config"
        self.name = "fixture_anaLambda"
        self.generator = self.root / "script" / GENERATOR.name
        self.generator.parent.mkdir(parents=True)
        shutil.copy2(str(GENERATOR), str(self.generator))
        self.analysis = self.config / "analysis" / "analysis_info_fixture.yaml"
        self.main = self.config / "mainconf" / ("main_" + self.name + ".yaml")
        self.write(self.analysis, (
            "analysis:\n  anaName: \"fixture_anaLambda\"\n"
            "  baseAnaMacro: \"anaGuardFixture\"\n"
            "  baseRunMacro: \"run_anaGuardFixture\"\n"
            "  mode: \"fxtmult\"\n"))
        if main_text is not None:
            self.write(self.main, main_text)
        templates = {
            "cuts/event/event.yaml": "minVz: -100\nmaxVz: 100\nmaxNTr: 0\n",
            "cuts/track/track.yaml": "minPt: 0.2\n",
            "cuts/pid/pid.yaml": "nSigmaProton: 2\n",
            "cuts/v0reco/v0.yaml": "maxDaughterDCA: 0.5\n",
            "cuts/mixing/mixing.yaml": "bufferSize: 20\n",
            "cuts/centrality/centrality.yaml": "enabled: true\nmode: refmult\n",
            "maker/maker_anaLambda.yaml": "nSigmaProton: 3\n",
            "hist/hist_anaLambda.yaml": "# fixture histogram\n",
            "mainconf/mainconf.yaml": (
                "event: cuts/event/event.yaml\ntrack: cuts/track/track.yaml\n"
                "pid: cuts/pid/pid.yaml\nv0: cuts/v0reco/v0.yaml\n"
                "mixing: cuts/mixing/mixing.yaml\ncentrality: cuts/centrality/centrality.yaml\n"
                "lambda: maker/maker___ANANAME__.yaml\n"
                "hist: hist/hist___ANANAME__.yaml\nanalysis: placeholder.yaml\n"),
        }
        for name, content in templates.items():
            self.write(self.config / name, content)
        # Existing concerns detect unauthorized mode updates, even without force.
        self.write(self.config / "cuts/event" / ("event_" + self.name + ".yaml"),
                   "minVz: -17\nmaxVz: 19\nmaxNTr: 7\n")
        self.write(self.config / "cuts/centrality" / ("centrality_" + self.name + ".yaml"),
                   "enabled: true\nmode: refmult\n")
        self.write(self.root / "analysis/anaPhi.C", "void anaPhi() {}\n")
        self.write(self.root / "analysis/run_anaPhi.C", "void run_anaPhi() {}\n")
        self.write(self.root / "analysis/anaGuardFixture.C", "// preserve existing analysis macro\n")
        self.write(self.root / "analysis/run_anaGuardFixture.C", "// preserve existing runner\n")

    @staticmethod
    def write(path, content):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content.encode("utf-8"))

    def snapshot(self):
        result = {}
        for path in self.root.rglob("*"):
            relative = str(path.relative_to(self.root))
            if path.is_file():
                result[relative] = (path.read_bytes(), path.stat().st_mtime_ns)
            elif path.is_dir():
                result[relative] = ("directory",)
        return result

    def run(self, force=False):
        command = [sys.executable, "-B", str(self.generator), str(self.analysis),
                   "--config-dir", str(self.config)]
        if force:
            command.append("--force")
        return subprocess.run(command, cwd=str(self.root), stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, universal_newlines=True)


class ManualConfigGuardTests(unittest.TestCase):
    def test_marked_cli_blocks_all_writes_with_and_without_force(self):
        for force in (False, True):
            with self.subTest(force=force), tempfile.TemporaryDirectory(prefix="star-manual-config-test-") as root:
                fixture = Fixture(root, MARKER + "\n# retain manually selected concerns\n")
                before = fixture.snapshot()
                result = fixture.run(force)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("manually maintained", result.stderr)
                self.assertIn("--force does not override", result.stderr)
                self.assertNotIn("Created:", result.stdout)
                self.assertEqual(before, fixture.snapshot())

    def test_marked_direct_mainconf_writer_is_also_guarded(self):
        with tempfile.TemporaryDirectory(prefix="star-manual-config-test-") as root:
            fixture = Fixture(root, MARKER + "\n")
            module = runpy.run_path(str(fixture.generator))
            before = fixture.snapshot()
            for force in (False, True):
                with self.subTest(force=force):
                    with self.assertRaisesRegex(RuntimeError, "manually maintained"):
                        module["write_mainconf"](str(fixture.config), fixture.name,
                                                 "analysis/analysis_info_fixture.yaml", force=force)
                    self.assertEqual(before, fixture.snapshot())

    def test_marker_can_be_after_other_lines_and_use_crlf(self):
        with tempfile.TemporaryDirectory(prefix="star-manual-config-test-") as root:
            fixture = Fixture(root, "# heading\r\n" + MARKER + "\r\n")
            before = fixture.snapshot()
            self.assertNotEqual(fixture.run().returncode, 0)
            self.assertEqual(before, fixture.snapshot())

    def test_unmarked_existing_mainconf_keeps_original_skip_force_behavior(self):
        for force in (False, True):
            with self.subTest(force=force), tempfile.TemporaryDirectory(prefix="star-manual-config-test-") as root:
                original = "# unmarked mainconf sentinel\n"
                fixture = Fixture(root, original)
                result = fixture.run(force)
                self.assertEqual(result.returncode, 0, result.stderr)
                if force:
                    self.assertIn("lambda:", fixture.main.read_text())
                    self.assertNotEqual(fixture.main.read_text(), original)
                else:
                    self.assertEqual(fixture.main.read_text(), original)
                self.assertTrue((fixture.config / "cuts/pid/pid_fixture_anaLambda.yaml").is_file())
                self.assertIn("mode: fxtmult",
                              (fixture.config / "cuts/centrality/centrality_fixture_anaLambda.yaml").read_text())
                self.assertIn("void anaGuardFixture",
                              (fixture.root / "analysis/anaGuardFixture.C").read_text())

    def test_similar_comments_are_not_the_exact_marker(self):
        comments = [
            "# config-generator: manual-extra",
            "# note: # config-generator: manual",
            "  # config-generator: manual",
            "# config-generator: manual ",
        ]
        for comment in comments:
            with self.subTest(comment=comment), tempfile.TemporaryDirectory(prefix="star-manual-config-test-") as root:
                fixture = Fixture(root, comment + "\n")
                result = fixture.run(force=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn("lambda:", fixture.main.read_text())

    def test_absent_mainconf_is_generated_as_before(self):
        with tempfile.TemporaryDirectory(prefix="star-manual-config-test-") as root:
            fixture = Fixture(root, None)
            result = fixture.run()
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue(fixture.main.is_file())
            self.assertIn("lambda:", fixture.main.read_text())


if __name__ == "__main__":
    unittest.main()
