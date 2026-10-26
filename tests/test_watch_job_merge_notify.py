"""Offline safety tests: no ROOT, Condor, email, or production file writes."""
import importlib.util
import json
from pathlib import Path
import subprocess
import tempfile
import types
import unittest
from unittest import mock
import xml.etree.ElementTree as ET

spec = importlib.util.spec_from_file_location(
    "watch_job_merge_notify", Path(__file__).resolve().parents[1] / "script/watch_job_merge_notify.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def xml_for(directory):
    return ET.fromstring('<job><output fromScratch="ana/ana_$JOBID.root" toURL="file:{}"/>'
                         '<stdout URL="file:{}/stdout.$JOBID.out"/></job>'.format(directory, directory))


def summary(missing=(), bad=(), incomplete=()):
    return dict(missingRoots=list(missing), badRootHeaders=list(bad),
                analysisCompletionUnconfirmed=list(incomplete))


class MergeNotifyTests(unittest.TestCase):
    def test_archive_failure_keeps_valid_merge_and_warns(self):
        with tempfile.TemporaryDirectory() as temp:
            state = Path(temp)
            result = dict(phase="merged", expectedSubjobs=1, checkedAt="test",
                          stateDirectory=temp, merge=dict(path="valid.root"))
            with mock.patch.object(module.archive_helper, "archive_merged_subjobs",
                                   side_effect=OSError("archive unavailable")):
                module.archive_after_merge(result, state / "runmeta.json", state)
            self.assertEqual(result["phase"], "merged")
            self.assertEqual(result["archive"]["phase"], "failed")
            msg = module.make_message("recipient@example.org", dict(anaName="test", jobid="A"*32), result)
            self.assertIn("要確認", str(msg["Subject"]))
            self.assertIn("valid.root", msg.get_content())
            self.assertIn("archive unavailable", msg.get_content())

    def test_failed_merge_is_never_archived(self):
        with mock.patch.object(module.archive_helper, "archive_merged_subjobs") as mover:
            module.archive_after_merge(dict(phase="failed"), Path("unused"), Path("unused"))
            mover.assert_not_called()

    def test_partial_requires_explicit_opt_in(self):
        with self.assertRaisesRegex(ValueError, "allow-partial"):
            module.select_inputs(xml_for("/tmp/test"), "ABC", range(3), summary(missing=[1]),
                                 {"exit_0": 2, "exit_84": 1}, False)
        paths, selected, excluded, partial, target = module.select_inputs(
            xml_for("/tmp/test"), "ABC", range(3), summary(missing=[1]),
            {"exit_0": 2, "exit_84": 1}, True)
        self.assertEqual(selected, [0, 2])
        self.assertEqual(excluded, [1])
        self.assertTrue(partial)
        self.assertEqual([p.name for p in paths], ["ana_ABC_0.root", "ana_ABC_2.root"])
        self.assertEqual(target.name, "ana_ABC_merge.root")

    def test_unconfirmed_analysis_is_excluded_even_when_root_exists(self):
        _, selected, excluded, partial, _ = module.select_inputs(
            xml_for("/tmp/test"), "ABC", range(3), summary(incomplete=[1]), {"exit_0": 3}, True)
        self.assertEqual((selected, excluded, partial), ([0, 2], [1], True))

    def test_bad_header_never_skipped_even_in_partial_mode(self):
        for allow in (False, True):
            with self.assertRaisesRegex(ValueError, "corrupt"):
                module.select_inputs(xml_for("/tmp/test"), "ABC", range(3), summary(bad=[1]),
                                     {"exit_0": 3}, allow)

    def test_empty_selection_is_error(self):
        with self.assertRaisesRegex(ValueError, "No completed"):
            module.select_inputs(xml_for("/tmp/test"), "ABC", range(1), summary(missing=[0]),
                                 {"exit_0": 1}, True)

    def test_numeric_order_is_retained_not_lexicographic(self):
        paths, selected, _, partial, _ = module.select_inputs(
            xml_for("/tmp/test"), "ABC", range(12), summary(), {"exit_0": 12}, False)
        self.assertEqual(selected, list(range(12)))
        self.assertEqual([p.name for p in paths][9:11], ["ana_ABC_9.root", "ana_ABC_10.root"])
        self.assertFalse(partial)

    def test_only_absolute_local_output_layout_is_allowed(self):
        for xml in [xml_for("relative/path"), ET.fromstring(
                '<job><output fromScratch="ana_$JOBID.root" toURL="root://example/root"/></job>')]:
            with self.assertRaises(ValueError):
                module.root_layout(xml, "ABC")

    def test_output_alias_is_resolved_for_container_inputs_and_destination(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            canonical = root / "gpfs-output"
            canonical.mkdir()
            alias = root / "star-u-output"
            alias.symlink_to(canonical, target_is_directory=True)
            directory, _, destination = module.root_layout(xml_for(alias), "ABC")
            self.assertEqual(directory, canonical)
            self.assertEqual(destination.parent, canonical)
            paths, _, _, _, target = module.select_inputs(
                xml_for(alias), "ABC", range(2), summary(), {"exit_0": 2}, False)
            self.assertEqual([p.parent for p in paths], [canonical, canonical])
            self.assertEqual(target, canonical / "ana_ABC_merge.root")

    def test_condor_completion_requires_all_events_and_empty_queue(self):
        self.assertFalse(module.notify.ready(2, {"1.0", "1.1"}, {"1.0": "exit_0"}, {}))
        self.assertFalse(module.notify.ready(1, {"1.0"}, {"1.0": "exit_0"}, {"1.0": 5}))
        self.assertTrue(module.notify.ready(1, {"1.0"}, {"1.0": "exit_84"}, {}))

    def test_hadd_nonzero_or_error_log_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            with mock.patch.object(module.subprocess, "run", return_value=types.SimpleNamespace(returncode=1)):
                with self.assertRaisesRegex(RuntimeError, "Command failed"):
                    module.run_logged(["hadd", "out.root"], directory / "rc.log")
            def root_error(command, stdout, stderr):
                stdout.write("Error in <TFile::ReadBuffer>: corrupt input\n")
                return types.SimpleNamespace(returncode=0)
            with mock.patch.object(module.subprocess, "run", side_effect=root_error):
                with self.assertRaisesRegex(RuntimeError, "ROOT error"):
                    module.run_logged(["hadd", "out.root"], directory / "error.log")

    def test_atomic_publish_refuses_existing_output(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source, target = root / "staging.root", root / "merged.root"
            source.write_bytes(b"root-new")
            target.write_bytes(b"root-existing")
            with self.assertRaises(FileExistsError):
                module.publish_no_clobber(source, target)
            self.assertEqual(target.read_bytes(), b"root-existing")
            other = root / "fresh.root"
            module.publish_no_clobber(source, other)
            self.assertEqual(other.read_bytes(), source.read_bytes())
            self.assertEqual(other.stat().st_ino, source.stat().st_ino)

    def test_histogram_validation_fails_on_missing_schema_or_sum(self):
        expected = {"h;1": ("TH1F:1:10:1:1", (3., 3., 9.))}
        module.compare_histograms(expected, expected)
        for actual in [{}, {"h;1": ("TH1D:1:10:1:1", (3., 3., 9.))},
                       {"h;1": ("TH1F:1:10:1:1", (3., 2., 9.))}]:
            with self.assertRaises(RuntimeError):
                module.compare_histograms(expected, actual)

    def test_audit_compiles_root5_macro_and_checks_report(self):
        with tempfile.TemporaryDirectory() as temp:
            state = Path(temp)
            inputs = [state / "a.root", state / "b.root"]
            def report(command, log):
                self.assertIn("merge_notify_audit.C+(", command[-1])
                (state / "input.tsv").write_text(
                    "".join("FILE\t{}\t1\n".format(p) for p in inputs)
                    + "HIST\th;1\tTH1F:1:10:1:1\t3\t3\t9\nDONE\t2\n")
            with mock.patch.object(module, "run_logged", side_effect=report):
                result = module.audit(inputs, "input", state, ["runner"])
            self.assertEqual(result["h;1"][1], (3., 3., 9.))
            with mock.patch.object(module, "run_logged"):
                (state / "invalid.tsv").write_text("DONE\t2\n")
                with self.assertRaisesRegex(RuntimeError, "input list mismatch"):
                    module.audit(inputs, "invalid", state, ["runner"])

    def test_merge_writes_manifest_batches_and_retains_originals(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            state = root / "state"
            state.mkdir()
            inputs = [root / "input_{}.root".format(i) for i in range(5)]
            for path in inputs:
                path.write_bytes(b"root-input")
            target = root / "merge.root"
            histogram = {"h;1": ("TH1F", (5., 5., 5.))}
            commands = []
            def hadd(command, log):
                commands.append(command)
                self.assertEqual(command[:2], ["runner", "hadd"])
                self.assertNotIn("-f", command)
                self.assertNotIn("-k", command)
                Path(command[2]).write_bytes(b"root-merged")
            with mock.patch.object(module, "audit", return_value=histogram), \
                    mock.patch.object(module, "run_logged", side_effect=hadd):
                result = module.merge_once(inputs, target, state, ["runner"], 2)
            self.assertEqual(result["inputCount"], 5)
            self.assertEqual(len(commands), 6)  # 5 -> 3 -> 2 -> 1
            self.assertEqual(target.read_bytes(), b"root-merged")
            self.assertTrue(all(p.read_bytes() == b"root-input" for p in inputs))
            manifest = json.loads((state / "input_manifest.json").read_text())
            self.assertEqual([entry["path"] for entry in manifest["inputs"]], [str(p) for p in inputs])
            self.assertEqual((state / "inputs.list").read_text().splitlines(), [str(p) for p in inputs])
            with self.assertRaisesRegex(RuntimeError, "already attempted"):
                module.merge_once(inputs, target, state, ["runner"], 2)

    def test_failed_audit_leaves_attempt_marker_and_no_published_output(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            state = root / "state"
            state.mkdir()
            source, target = root / "source.root", root / "target.root"
            source.write_bytes(b"root-input")
            with mock.patch.object(module, "audit", side_effect=RuntimeError("corrupt ROOT")), \
                    mock.patch.object(module, "run_logged") as hadd:
                with self.assertRaisesRegex(RuntimeError, "corrupt ROOT"):
                    module.merge_once([source], target, state, ["runner"], 2)
            self.assertTrue((state / "merge_attempt.json").exists())
            self.assertFalse(target.exists())
            hadd.assert_not_called()

    def test_mutated_input_is_not_published(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            state = root / "state"
            state.mkdir()
            source, target = root / "source.root", root / "target.root"
            source.write_bytes(b"root-input")
            histogram = {"h;1": ("TH1F", (1., 1., 1.))}
            def hadd(command, log):
                Path(command[2]).write_bytes(b"root-merged")
                source.write_bytes(b"root-input-changed")
            with mock.patch.object(module, "audit", return_value=histogram), \
                    mock.patch.object(module, "run_logged", side_effect=hadd):
                with self.assertRaisesRegex(RuntimeError, "changed"):
                    module.merge_once([source], target, state, ["runner"], 2)
            self.assertFalse(target.exists())

    def test_failure_or_partial_merge_is_warning_email(self):
        meta = dict(anaName="test", jobid="A" * 32)
        base = dict(expectedSubjobs=3, checkedAt="test", stateDirectory="/tmp/test")
        for phase, partial, warning in [("merged", False, False), ("merged", True, True),
                                        ("failed", False, True)]:
            result = dict(base, phase=phase, partial=partial)
            msg = module.make_message("recipient@example.org", meta, result)
            self.assertEqual("要確認" in str(msg["Subject"]), warning)
            self.assertEqual(msg["To"], "recipient@example.org")

    def test_finish_without_email_never_invokes_sendmail(self):
        with tempfile.TemporaryDirectory() as temp, mock.patch.object(module.notify, "send_once") as sender:
            result = dict(phase="merged")
            self.assertEqual(module.finish(result, {}, types.SimpleNamespace(email=None), Path(temp)), 0)
            sender.assert_not_called()

    def test_merge_notification_does_not_retry_ambiguous_mail(self):
        with tempfile.TemporaryDirectory() as temp, mock.patch.object(module.notify.subprocess, "run") as sender:
            state = Path(temp)
            sender.side_effect = subprocess.TimeoutExpired("sendmail", 45)
            args = types.SimpleNamespace(email="recipient@example.org")
            meta = dict(anaName="test", jobid="A" * 32)
            result = dict(phase="failed", expectedSubjobs=3, checkedAt="test", stateDirectory=str(state))
            with self.assertRaises(SystemExit):
                module.finish(result, meta, args, state)
            self.assertTrue((state / "delivery_attempt.json").exists())
            with self.assertRaises(SystemExit):
                module.finish(result, meta, args, state)
            sender.assert_called_once()

    def test_strict_selection_failure_keeps_excluded_indices_in_warning(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            jobid = "A" * 32
            prefix = "test" + jobid
            run = root / "job/run"
            run.mkdir(parents=True)
            for index in range(3):
                (run / "{}_{}.list".format(prefix, index)).write_text("input.root\n")
            xml_path = root / "job.xml"
            xml_path.write_text(ET.tostring(xml_for(root), encoding="unicode"))
            mainconf, runner = root / "main.yaml", root / "runner.sh"
            mainconf.touch()
            runner.touch()
            events = root / "events.log"
            events.touch()
            meta = dict(anaName="test", jobid=jobid, jobPrefix=prefix, projectRoot=str(root),
                        submit=dict(host=module.socket.getfqdn(), mainconf=str(mainconf),
                                    joblistSnapshot=dict(path=str(xml_path))))
            metadata = root / "runmeta.json"
            metadata.write_text(json.dumps(meta))
            argv = ["watch_job_merge_notify.py", "--runmeta", str(metadata),
                    "--event-log", str(events), "--root-runner", str(runner),
                    "--no-email", "--watch", "--settle-sec", "0"]
            submitted = {"1.0", "1.1", "1.2"}
            with mock.patch.object(module.sys, "argv", argv), \
                    mock.patch.object(module.notify, "parse_events", return_value=(
                        submitted, {key: "exit_0" for key in submitted})), \
                    mock.patch.object(module.notify, "queue_states", return_value={}), \
                    mock.patch.object(module.notify, "output_summary", return_value=summary(
                        missing=[1], incomplete=[1, 2])), \
                    mock.patch.object(module, "merge_once") as merger, \
                    mock.patch.object(module, "finish", return_value=1) as finish:
                self.assertEqual(module.main(), 1)
            result = finish.call_args[0][0]
            self.assertEqual(result["phase"], "failed")
            self.assertEqual(result["excludedIndices"], [1, 2])
            msg = module.make_message("recipient@example.org", meta, result)
            self.assertIn("除外したsubjob番号: [1, 2]", msg.get_content())
            merger.assert_not_called()


if __name__ == "__main__":
    unittest.main()
