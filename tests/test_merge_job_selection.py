"""Explicit common-population merging without ROOT, Condor or email actions."""
import contextlib
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock
import xml.etree.ElementTree as ET

spec = importlib.util.spec_from_file_location(
    "merge_job_selection", Path(__file__).resolve().parents[1] / "script/watch_job_merge_notify.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def xml_for(directory):
    return ET.fromstring('<job><output fromScratch="ana/ana_$JOBID.root" toURL="file:{}"/>'
                         '<stdout URL="file:{}/stdout.$JOBID.out"/></job>'.format(directory, directory))


def summary(missing=(), bad=(), incomplete=()):
    return dict(missingRoots=list(missing), badRootHeaders=list(bad),
                analysisCompletionUnconfirmed=list(incomplete))


class JobSelectionTests(unittest.TestCase):
    def test_selection_records_original_bytes_and_sorts_numerically(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "common.list"
            raw = b"10\n2\n0\n"
            path.write_bytes(raw)
            selected = module.load_job_selection(path, range(11))
            self.assertEqual(selected, dict(path=str(path.resolve()),
                sha256=hashlib.sha256(raw).hexdigest(), jobIds=[0, 2, 10], count=3))

    def test_empty_duplicate_unknown_and_noncanonical_ids_are_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "common.list"
            for text in ["", "\n", "0\n\n1\n", " 1\n", "1 2\n", "1.0\n", "+1\n",
                         "-1\n", "01\n", "#1\n", "1\n1\n", "99\n"]:
                with self.subTest(text=text):
                    path.write_text(text)
                    with self.assertRaises(ValueError):
                        module.load_job_selection(path, range(3))

    def test_default_population_and_order_are_preserved(self):
        paths, selected, excluded, partial, _ = module.select_inputs(
            xml_for("/tmp/test-selection"), "ABC", range(12), summary(), {"exit_0": 12}, False)
        self.assertEqual(selected, list(range(12)))
        self.assertEqual(excluded, [])
        self.assertFalse(partial)
        self.assertEqual(paths[10].name, "ana_ABC_10.root")

    def test_valid_ep1_job_is_excluded_by_population_not_bad_root(self):
        missing = list(range(34))
        include = [i for i in range(4058) if i not in set(missing + [3496])]
        output = summary(missing=missing, incomplete=missing)
        paths, selected, excluded, partial, _ = module.select_inputs(
            xml_for("/tmp/test-selection"), "ABC", range(4058), output,
            {"exit_0": 4058}, True, include)
        reasons = module.exclusion_reasons(range(4058), output, include)
        self.assertEqual(len(paths), 4023)
        self.assertEqual(selected, include)
        self.assertEqual(excluded, missing + [3496])
        self.assertTrue(partial)
        self.assertEqual(reasons["3496"], ["outside_requested_selection"])
        self.assertEqual(reasons["0"], ["missing_root", "analysis_completion_unconfirmed"])

    def test_subset_still_requires_partial_opt_in_and_does_not_skip_corruption(self):
        with self.assertRaisesRegex(ValueError, "allow-partial"):
            module.select_inputs(xml_for("/tmp/test-selection"), "ABC", range(3),
                                 summary(), {"exit_0": 3}, False, [0, 2])
        with self.assertRaisesRegex(ValueError, "corrupt"):
            module.select_inputs(xml_for("/tmp/test-selection"), "ABC", range(3),
                                 summary(bad=[1]), {"exit_0": 3}, True, [0, 2])
        with self.assertRaisesRegex(ValueError, "No completed"):
            module.select_inputs(xml_for("/tmp/test-selection"), "ABC", range(2),
                                 summary(missing=[0]), {"exit_0": 2}, True, [0])

    def test_modified_selection_fails_before_a_merge_attempt(self):
        with tempfile.TemporaryDirectory() as temp:
            state = Path(temp)
            include = state / "common.list"
            include.write_text("0\n")
            selection = module.load_job_selection(include, range(2))
            include.write_text("1\n")
            with self.assertRaisesRegex(ValueError, "changed"):
                module.merge_once([], state / "output.root", state, ["runner"], 2,
                                  selection=selection)
            self.assertFalse((state / "merge_attempt.json").exists())

    def test_selected_population_is_saved_with_input_provenance(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            state = root / "state"
            state.mkdir()
            include = root / "common.list"
            include.write_text("2\n0\n")
            selection = module.load_job_selection(include, range(3))
            inputs = [root / "input_0.root", root / "input_2.root"]
            for path in inputs:
                path.write_bytes(b"root-input")
            histogram = {"h;1": ("TH1F", (2., 2., 2.))}
            def hadd(command, log):
                Path(command[2]).write_bytes(b"root-merged")
            with mock.patch.object(module, "audit", return_value=histogram), \
                    mock.patch.object(module, "run_logged", side_effect=hadd):
                module.merge_once(inputs, root / "merge.root", state, ["runner"], 2,
                                  selection=selection)
            for name in ["input_manifest.json", "merge_attempt.json"]:
                record = json.loads((state / name).read_text())
                self.assertEqual(record["includeJobIds"], selection)
                self.assertEqual(record["selectedInputCount"], 2)
            self.assertEqual((state / "inputs.list").read_text().splitlines(),
                             [str(path) for path in inputs])

    def fixture(self, root, watch=True):
        jobid = "A" * 32
        prefix = "test" + jobid
        run = root / "job/run"
        run.mkdir(parents=True)
        for index in range(4):
            (run / "{}_{}.list".format(prefix, index)).write_text("input.root\n")
        xml_path = root / "job.xml"
        xml_path.write_text(ET.tostring(xml_for(root), encoding="unicode"))
        mainconf, runner, events = root / "main.yaml", root / "runner.sh", root / "events.log"
        for path in [mainconf, runner, events]:
            path.touch()
        meta = dict(anaName="test", jobid=jobid, jobPrefix=prefix, projectRoot=str(root),
                    submit=dict(host=module.socket.getfqdn(), mainconf=str(mainconf),
                                joblistSnapshot=dict(path=str(xml_path))))
        metadata = root / "runmeta.json"
        metadata.write_text(json.dumps(meta))
        include = root / "common.list"
        include.write_text("2\n0\n")
        argv = ["watch_job_merge_notify.py", "--runmeta", str(metadata),
                "--event-log", str(events), "--root-runner", str(runner), "--no-email",
                "--watch" if watch else "--check", "--settle-sec", "0", "--allow-partial",
                "--include-job-ids", str(include)]
        submitted = {"1.{}".format(i) for i in range(4)}
        return argv, include, metadata, xml_path, run / "merge_notify" / ("test_" + jobid), submitted

    def test_check_only_reports_the_frozen_population_without_creating_state(self):
        with tempfile.TemporaryDirectory() as temp:
            argv, include, _, _, state, submitted = self.fixture(Path(temp), watch=False)
            output = io.StringIO()
            with mock.patch.object(module.sys, "argv", argv), \
                    mock.patch.object(module.notify, "parse_events", return_value=(
                        submitted, {key: "exit_0" for key in submitted})), \
                    mock.patch.object(module.notify, "queue_states", return_value={}), \
                    contextlib.redirect_stdout(output):
                self.assertEqual(module.main(), 0)
            record = json.loads(output.getvalue())
            self.assertEqual(record["includeJobIds"]["jobIds"], [0, 2])
            self.assertEqual(record["expectedSubjobs"], 4)
            self.assertFalse(state.exists())

    def test_cli_passes_selection_and_records_the_same_population_everywhere(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            argv, include, metadata, xml_path, state, submitted = self.fixture(root)
            original_meta, original_xml = metadata.read_bytes(), xml_path.read_bytes()
            with mock.patch.object(module.sys, "argv", argv), \
                    mock.patch.object(module.notify, "parse_events", return_value=(
                        submitted, {key: "exit_0" for key in submitted})), \
                    mock.patch.object(module.notify, "queue_states", return_value={}), \
                    mock.patch.object(module.notify, "output_summary", return_value=summary(
                        missing=[1], incomplete=[1])), \
                    mock.patch.object(module, "merge_once", return_value=dict(path="valid.root")) as merger, \
                    mock.patch.object(module, "archive_after_merge"), \
                    mock.patch.object(module, "finish", return_value=0), \
                    contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(module.main(), 0)
            for name in ["settings.json", "result.json", "status.json"]:
                record = json.loads((state / name).read_text())
                self.assertEqual(record["includeJobIds"]["jobIds"], [0, 2])
                self.assertEqual(record["includeJobIds"]["count"], 2)
            result = json.loads((state / "result.json").read_text())
            self.assertEqual(result["selectedIndices"], [0, 2])
            self.assertEqual(result["exclusionReasons"]["3"], ["outside_requested_selection"])
            self.assertEqual(merger.call_args.kwargs["selection"], result["includeJobIds"])
            self.assertEqual(metadata.read_bytes(), original_meta)
            self.assertEqual(xml_path.read_bytes(), original_xml)


if __name__ == "__main__":
    unittest.main()
