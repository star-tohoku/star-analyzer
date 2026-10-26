import importlib.util
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


class NotifyTests(unittest.TestCase):
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
