"""Offline archive safety tests. Only temporary test ROOT-like files are moved."""
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest import mock

spec = importlib.util.spec_from_file_location(
    "archive_merged_subjobs", Path(__file__).resolve().parents[1] / "script/archive_merged_subjobs.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class ArchiveTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.jobid, self.ana = "A" * 32, "test"
        self.out = self.root / "rootfile/test"
        self.out.mkdir(parents=True)
        self.state = self.root / "job/run/merge_notify" / (self.ana + "_" + self.jobid)
        self.state.mkdir(parents=True)
        self.inputs = [self.out / ("test_" + self.jobid + "_" + str(i) + ".root") for i in (0, 2)]
        for path in self.inputs:
            path.write_bytes(b"root-input")
        self.failed = self.out / ("test_" + self.jobid + "_1.root")
        self.failed.write_bytes(b"failed-unmerged-data")
        self.staged = self.out / ".staging/round.root"
        self.staged.parent.mkdir()
        self.staged.write_bytes(b"root-merged")
        self.merged = self.out / ("test_" + self.jobid + "_merge.root")
        os.link(str(self.staged), str(self.merged))
        self.xml = self.root / "job/job.xml"
        self.xml.write_text('<job><output fromScratch="test_$JOBID.root" toURL="file:{}"/></job>'.format(self.out))
        self.meta = self.root / "runmeta.json"
        self.write(self.meta, dict(anaName=self.ana, jobid=self.jobid, projectRoot=str(self.root),
                                  submit=dict(joblistSnapshot=dict(path=str(self.xml)))))
        self.result = dict(phase="merged", stateDirectory=str(self.state), selectedIndices=[0, 2], expectedSubjobs=3,
                           excludedIndices=[1], partial=True, merge=dict(
                               path=str(self.merged), size=self.merged.stat().st_size,
                               histogramCount=1, validation="ROOT recursive read + additive histogram fingerprints",
                               inputCount=2, publishedAt=module.now()))
        self.write(self.state / "result.json", self.result)
        self.manifest = dict(inputs=[dict(path=str(p), **module.regular(p)) for p in self.inputs])
        self.write(self.state / "input_manifest.json", self.manifest)
        (self.state / "merged_audit.tsv").write_text("FILE\t{}\t1\nHIST\th;1\tTH1D\t2\t2\t2\nDONE\t1\n".format(self.staged))
        self.target = self.out / "subjobs" / self.jobid

    def write(self, path, data):
        path.write_text(json.dumps(data))

    def archive(self):
        return module.archive_merged_subjobs(self.meta)

    def assert_sources_retained(self):
        self.assertTrue(all(p.exists() for p in self.inputs))

    def test_success_partial_only_and_historical_records_unchanged(self):
        original = {p: p.read_bytes() for p in (self.meta, self.state / "result.json", self.state / "input_manifest.json")}
        result = self.archive()
        self.assertEqual(result["phase"], "archived")
        for path, entry in zip(self.inputs, self.manifest["inputs"]):
            self.assertFalse(path.exists())
            module.verify_fingerprint(self.target / path.name, entry)
        self.assertEqual(self.failed.read_bytes(), b"failed-unmerged-data")
        self.assertEqual(self.merged.read_bytes(), b"root-merged")
        for path, content in original.items():
            self.assertEqual(path.read_bytes(), content)
        self.assertTrue((self.target / "README.md").is_file())
        self.assertEqual(self.archive()["phase"], "archived")

    def test_resume_after_link_before_unlink(self):
        unlink = module.os.unlink
        def interrupt(path, *args, **kwargs):
            if path == str(self.inputs[0]):
                raise OSError("simulated interruption")
            return unlink(path, *args, **kwargs)
        with mock.patch.object(module.os, "unlink", side_effect=interrupt):
            with self.assertRaisesRegex(OSError, "interruption"):
                self.archive()
        self.assertTrue(self.inputs[0].exists())
        self.assertTrue((self.target / self.inputs[0].name).exists())
        self.assertEqual(self.archive()["phase"], "archived")
        self.assertFalse(self.inputs[0].exists())

    def test_resume_after_original_unlink(self):
        unlink = module.os.unlink
        def interrupt(path, *args, **kwargs):
            result = unlink(path, *args, **kwargs)
            if path == str(self.inputs[0]):
                raise OSError("interruption after unlink")
            return result
        with mock.patch.object(module.os, "unlink", side_effect=interrupt):
            with self.assertRaisesRegex(OSError, "after unlink"):
                self.archive()
        self.assertFalse(self.inputs[0].exists())
        self.assertEqual((self.target / self.inputs[0].name).read_bytes(), b"root-input")
        self.assertEqual(self.archive()["phase"], "archived")

    def test_collision_prevalidated_before_any_input_moves(self):
        self.target.mkdir(parents=True)
        (self.target / self.inputs[1].name).write_bytes(b"do-not-overwrite")
        with self.assertRaisesRegex(FileExistsError, "collision"):
            self.archive()
        self.assert_sources_retained()

    def test_mutated_input_prevalidated_before_moves(self):
        self.inputs[1].write_bytes(b"mutated")
        with self.assertRaisesRegex(ValueError, "fingerprint"):
            self.archive()
        self.assert_sources_retained()

    def test_missing_merge(self):
        self.merged.unlink()
        with self.assertRaises(FileNotFoundError):
            self.archive()
        self.assert_sources_retained()

    def test_missing_input_aborts_before_other_inputs_move(self):
        self.inputs[1].unlink()
        with self.assertRaisesRegex(FileNotFoundError, "Missing merge input"):
            self.archive()
        self.assertTrue(self.inputs[0].is_file())

    def test_mutated_archive_is_not_treated_as_completed_resume(self):
        self.archive()
        (self.target / self.inputs[0].name).write_bytes(b"changed-archive")
        with self.assertRaisesRegex(ValueError, "fingerprint"):
            self.archive()

    def test_merge_modified_in_place_after_publication(self):
        self.merged.write_bytes(b"root-altered")
        later = module.datetime.datetime.fromisoformat(self.result["merge"]["publishedAt"]).timestamp() + 2
        os.utime(str(self.merged), (later, later))
        self.result["merge"]["size"] = self.merged.stat().st_size
        self.write(self.state / "result.json", self.result)
        with self.assertRaisesRegex(ValueError, "after validation"):
            self.archive()
        self.assert_sources_retained()

    def test_publication_fingerprint_rejects_modified_merge(self):
        self.result["merge"]["fingerprint"] = module.regular(self.merged)
        self.write(self.state / "result.json", self.result)
        os.utime(str(self.merged), (1, 1))
        with self.assertRaisesRegex(ValueError, "fingerprint"):
            self.archive()
        self.assert_sources_retained()

    def test_corrupt_published_root_header(self):
        with self.merged.open("r+b") as stream:
            stream.write(b"xxxx")
        with self.assertRaisesRegex(ValueError, "ROOT header"):
            self.archive()
        self.assert_sources_retained()

    def test_wrong_job_scope(self):
        self.manifest["inputs"][1]["path"] = str(self.failed)
        self.write(self.state / "input_manifest.json", self.manifest)
        with self.assertRaisesRegex(ValueError, "scope"):
            self.archive()
        self.assert_sources_retained()

    def test_not_merged(self):
        self.result["phase"] = "failed"
        self.write(self.state / "result.json", self.result)
        with self.assertRaisesRegex(ValueError, "phase=merged"):
            self.archive()
        self.assert_sources_retained()

    def test_symlink_destination(self):
        elsewhere = self.root / "elsewhere"
        elsewhere.mkdir()
        (self.out / "subjobs").symlink_to(elsewhere)
        with self.assertRaisesRegex(ValueError, "Symlink"):
            self.archive()
        self.assert_sources_retained()

    def test_symlink_input(self):
        saved = self.root / "original.root"
        self.inputs[1].rename(saved)
        self.inputs[1].symlink_to(saved)
        with self.assertRaisesRegex(ValueError, "Symlink"):
            self.archive()
        self.assertTrue(self.inputs[0].exists())

    def test_unvalidated_merge_or_bad_count(self):
        self.result["merge"]["inputCount"] = 3
        self.write(self.state / "result.json", self.result)
        with self.assertRaisesRegex(ValueError, "count"):
            self.archive()
        self.assert_sources_retained()

    def test_replaced_merge_rejected_even_same_bytes(self):
        self.merged.unlink()
        self.merged.write_bytes(b"root-merged")
        with self.assertRaisesRegex(ValueError, "staging"):
            self.archive()
        self.assert_sources_retained()


if __name__ == "__main__":
    unittest.main()
