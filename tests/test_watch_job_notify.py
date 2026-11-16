import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock
import xml.etree.ElementTree as ET

spec = importlib.util.spec_from_file_location("watch_job_notify", Path(__file__).resolve().parents[1] / "script/watch_job_notify.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def event(code, job="10.000.000", body=""):
    return "{} ({}) 2026-10-01 01:00:00 event\n{}\n...\n".format(code, job, body)


def message_fixture(good=True):
    meta = {"anaName": "PRIVATE_ANALYSIS_SENTINEL", "jobid": "0123456789ABCDEF0123456789ABCDEF"}
    status = {"expectedSubjobs": 26994, "outcomes": {"exit_0": 26994},
              "checkedAt": "PRIVATE_CHECK_TIME", "stateDirectory": "/private/internal-monitor",
              "output": {"rootsWithValidHeader": 26994, "rootDirectory": "/private/internal-output",
                         "missingRoots": [], "badRootHeaders": [],
                         "analysisCompletionUnconfirmed": [] if good else [17]}}
    return meta, status


class NotifyTests(unittest.TestCase):
    def assert_minimal_message(self, message, meta):
        visible = str(message["Subject"]) + "\n" + message.get_content()
        for private in (meta["anaName"], meta["jobid"], "/private/", "26994", "17",
                        "SUMS", "quota", "PRIVATE_CHECK_TIME"):
            self.assertNotIn(private, visible)
        self.assertNotIn(meta["jobid"], str(message["Message-ID"]))
        self.assertEqual(str(message["To"]), "recipient@example.org")

    def test_minimal_success_has_no_run_diagnostics(self):
        meta, status = message_fixture()
        with mock.patch.object(module.email.utils, "make_msgid", return_value="<opaque@example.org>"):
            message = module.make_message("recipient@example.org", meta, status, content="minimal")
        self.assert_minimal_message(message, meta)
        self.assertIn("完了", str(message["Subject"]))
        self.assertNotIn("要確認", str(message["Subject"]))
        self.assertIn("出力回収と解析完了ログを確認", message.get_content())
        self.assertEqual(str(message["Message-ID"]), "<opaque@example.org>")

    def test_minimal_warning_has_no_run_diagnostics(self):
        meta, status = message_fixture(good=False)
        with mock.patch.object(module.email.utils, "make_msgid", return_value="<opaque@example.org>"):
            message = module.make_message("recipient@example.org", meta, status, content="minimal")
        self.assert_minimal_message(message, meta)
        self.assertIn("要確認", str(message["Subject"]))
        self.assertIn("正常完了を確認できない", message.get_content())
        self.assertIn("このチャット", message.get_content())

    def test_detailed_is_the_backward_compatible_default(self):
        meta, status = message_fixture()
        default = module.make_message("recipient@example.org", meta, status)
        explicit = module.make_message("recipient@example.org", meta, status, content="detailed")
        self.assertEqual(default.get_content(), explicit.get_content())
        self.assertEqual(default["Subject"], explicit["Subject"])
        self.assertEqual(default["Message-ID"], explicit["Message-ID"])
        self.assertIn(meta["anaName"], str(default["Subject"]))
        self.assertIn(meta["jobid"], default.get_content())
        self.assertIn(status["output"]["rootDirectory"], default.get_content())

    def test_unknown_email_content_is_rejected(self):
        meta, status = message_fixture()
        with self.assertRaises(ValueError):
            module.make_message("recipient@example.org", meta, status, content="unknown")

    def test_minimal_cli_records_content_and_sends_only_minimal_message(self):
        self.check_minimal_cli_state(custom_state=False)

    def test_custom_state_dir_is_read_only_for_check_and_used_for_watch(self):
        self.check_minimal_cli_state(custom_state=True)

    def check_minimal_cli_state(self, custom_state):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = root / "job/run"
            run.mkdir(parents=True)
            (run / "PRIVATE_PREFIX_0.list").write_text("input.root\n")
            snapshot = root / "job.xml"
            snapshot.write_text("<job/>")
            meta, fixture_status = message_fixture()
            meta.update(projectRoot=str(root), jobPrefix="PRIVATE_PREFIX",
                        submit={"host": "test-submit-host", "joblistSnapshot": {"path": str(snapshot)}})
            runmeta = root / "runmeta.json"
            runmeta.write_text(json.dumps(meta))
            event_log = root / "events.log"
            event_log.write_text(event("000") + event("005", body="Normal termination (return value 0)"))
            output = fixture_status["output"]
            output["rootsWithValidHeader"] = 1
            argv = ["watch_job_notify.py", "--runmeta", str(runmeta), "--event-log", str(event_log),
                    "--email", "recipient@example.org", "--email-content", "minimal", "--watch", "--settle-sec", "0"]
            default_state = run / "notify" / (meta["anaName"] + "_" + meta["jobid"])
            state = root / "custom-monitor-state" if custom_state else default_state
            if custom_state:
                argv.extend(["--state-dir", str(state / ".." / state.name)])
            with mock.patch("sys.argv", argv), \
                    mock.patch.object(module.socket, "getfqdn", return_value="test-submit-host"), \
                    mock.patch.object(module.fcntl, "flock", side_effect=lambda handle, operation: self.addCleanup(handle.close)), \
                    mock.patch.object(module, "queue_states", return_value={}), \
                    mock.patch.object(module, "output_summary", return_value=output), \
                    mock.patch.object(module.subprocess, "run") as sender, \
                    mock.patch("builtins.print"):
                if custom_state:
                    argv[argv.index("--watch")] = "--check"
                    module.main()
                    self.assertFalse(state.exists())
                    sender.assert_not_called()
                    argv[argv.index("--check")] = "--watch"
                module.main()
            for name in ("settings.json", "status.json", "sent.json"):
                record = json.loads((state / name).read_text())
                self.assertEqual(record["emailContent"], "minimal")
                self.assertEqual(record["stateDirectory"], str(state.resolve()))
            self.assertTrue((state / "message.eml").exists())
            if custom_state:
                self.assertFalse(default_state.exists())
            sender.assert_called_once()
            sent = module.email.message_from_bytes(sender.call_args.kwargs["input"], policy=module.email.policy.default)
            self.assert_minimal_message(sent, meta)
            self.assertIn("出力回収と解析完了ログを確認", sent.get_content())

    def test_send_once_success_prevents_duplicate(self):
        with tempfile.TemporaryDirectory() as directory, mock.patch.object(module.subprocess, "run") as sender:
            state = Path(directory)
            message = module.email.message.EmailMessage()
            message["Message-ID"] = "<test@example.org>"
            module.send_once(state, message, {}, "recipient@example.org")
            self.assertTrue((state / "sent.json").exists())
            with self.assertRaises(SystemExit):
                module.send_once(state, message, {}, "recipient@example.org")
            sender.assert_called_once()

    def test_send_timeout_blocks_automatic_retry(self):
        with tempfile.TemporaryDirectory() as directory, mock.patch.object(module.subprocess, "run") as sender:
            sender.side_effect = module.subprocess.TimeoutExpired("sendmail", 45)
            state = Path(directory)
            message = module.email.message.EmailMessage()
            message["Message-ID"] = "<test@example.org>"
            with self.assertRaises(SystemExit):
                module.send_once(state, message, {}, "recipient@example.org")
            self.assertTrue((state / "delivery_attempt.json").exists())
            self.assertFalse((state / "sent.json").exists())
            with self.assertRaises(SystemExit):
                module.send_once(state, message, {}, "recipient@example.org")
            sender.assert_called_once()

    def test_record_failure_after_send_stops_instead_of_retrying(self):
        with tempfile.TemporaryDirectory() as directory, mock.patch.object(module.subprocess, "run") as sender:
            state = Path(directory)
            message = module.email.message.EmailMessage()
            message["Message-ID"] = "<test@example.org>"
            save = module.save_json
            def fail_sent_record(path, data):
                if path.name == "sent.json":
                    raise OSError("storage unavailable")
                save(path, data)
            with mock.patch.object(module, "save_json", side_effect=fail_sent_record):
                with self.assertRaises(SystemExit):
                    module.send_once(state, message, {}, "recipient@example.org")
            with self.assertRaises(SystemExit):
                module.send_once(state, message, {}, "recipient@example.org")
            sender.assert_called_once()

    def test_complete_and_trailing_resource_event(self):
        text = event("000") + event("005", body="Normal termination (return value 0)") + event("006")
        submitted, terminal = module.parse_events(text)
        self.assertEqual(terminal, {"10.0": "exit_0"})
        self.assertTrue(module.ready(1, submitted, terminal, {}))

    def test_partial_terminal_is_not_complete(self):
        text = event("000") + "005 (10.000.000) 2026-10-01 event\nNormal termination (return value 0)"
        submitted, terminal = module.parse_events(text)
        self.assertFalse(module.ready(1, submitted, terminal, {}))

    def test_failure_removed_and_signal(self):
        for code, body, expected in [("005", "Normal termination (return value 84)", "exit_84"),
                                     ("005", "Abnormal termination (signal 11)", "signal_or_unknown"),
                                     ("009", "removed", "removed")]:
            submitted, terminal = module.parse_events(event("000") + event(code, body=body))
            self.assertTrue(module.ready(1, submitted, terminal, {}))
            self.assertEqual(terminal["10.0"], expected)

    def test_queue_pending_and_missing_expected(self):
        submitted, terminal = module.parse_events(event("000") + event("005", body="Normal termination (return value 0)"))
        self.assertFalse(module.ready(2, submitted, terminal, {}))
        self.assertFalse(module.ready(1, submitted, terminal, {"10.0": 5}))

    def test_reexecution_invalidates_terminal(self):
        text = event("000") + event("005", body="Normal termination (return value 0)") + event("001")
        submitted, terminal = module.parse_events(text)
        self.assertFalse(module.ready(1, submitted, terminal, {}))

    def test_outputs_do_not_count_unrelated_files(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "ana_ABC_0.root").write_bytes(b"root1234")
            (root / "ana_ABC_1.root").write_bytes(b"bad!1234")
            (root / "ana_ABC_merge.root").write_bytes(b"root1234")
            (root / "stdout.ABC_0.out").write_text("FemtoLambda requested=3 completed=3 reconstructed=2 status=0\n")
            xml = ET.fromstring('<job><output fromScratch="ana/ana_$JOBID.root" toURL="file:{}"/><stdout URL="file:{}/stdout.$JOBID.out"/></job>'.format(root, root))
            report = module.output_summary(xml, "ABC", [0, 1, 2])
            self.assertEqual(report["rootsWithValidHeader"], 1)
            self.assertEqual(report["badRootHeaders"], [1])
            self.assertEqual(report["missingRoots"], [2])
            self.assertEqual(report["analysisCompletionUnconfirmed"], [1, 2])

    def test_warning_email_when_analysis_unconfirmed(self):
        meta = {"anaName": "test", "jobid": "A" * 32}
        status = {"expectedSubjobs": 1, "outcomes": {"exit_0": 1}, "checkedAt": "test", "stateDirectory": "/tmp/test",
                  "output": {"rootsWithValidHeader": 1, "rootDirectory": "/tmp", "missingRoots": [],
                             "badRootHeaders": [], "analysisCompletionUnconfirmed": [0]}}
        msg = module.make_message("recipient@example.org", meta, status)
        self.assertIn("要確認", str(msg["Subject"]))
        self.assertIn("正常完了は確認できません", msg.get_content())
        self.assertEqual(str(msg["To"]), "recipient@example.org")


if __name__ == "__main__":
    unittest.main()
