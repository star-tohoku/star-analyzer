#!/usr/bin/env python3
"""Static entry/config contract; no STAR import or input-file access required."""
from pathlib import Path
import subprocess
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPECIES = {"d": "deuteron", "t": "triton", "3He": "he3", "4He": "he4"}


def flat(path):
    result = {}
    for line in path.read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        key, value = line.split(":", 1)
        if key in result:
            raise AssertionError("duplicate YAML key: " + key)
        result[key.strip()] = value.strip().strip('"\'')
    return result


class FemtoLambdaEntries(unittest.TestCase):
    def test_eight_mainconfs_reference_presets_and_real_macros(self):
        for energy in ("13p5", "3p85"):
            for short, species in SPECIES.items():
                name = f"auau{energy}_anaFemtoLambda_{short}_KFParticle_highpurity"
                main = flat(ROOT / f"config/mainconf/main_{name}.yaml")
                self.assertEqual(set(main), {"analysis", "event", "centrality", "kf", "nuclearid",
                                             "maker", "mixing", "femtoHist", "nuclearHist"})
                self.assertEqual(main["kf"], f"cuts/kf/kf_auau{energy}_anaLambda_KFParticle_highpurity.yaml")
                for path in main.values():
                    self.assertTrue((ROOT / "config" / path).is_file(), path)
                info = (ROOT / "config" / main["analysis"]).read_text()
                self.assertIn(f'baseRunMacro: "run_anaFemtoLambda_{short}"', info)
                self.assertIn(f'baseAnaMacro: "anaFemtoLambda_{short}"', info)
                self.assertIn('mode: "fxtmult"', info)
                self.assertIn("fixedTarget_" + ("2020" if energy == "13p5" else "2019"), info)
                maker = flat(ROOT / "config" / main["maker"])
                self.assertEqual(maker["speciesKeys"], "lambda," + species)
                self.assertEqual(maker["species_lambda_particleKey"], "lambda_kf")
                self.assertEqual(maker["lambdaEventPolicy"], "lambda_imp5_compat")
                self.assertEqual(maker["nuclearSelectionProfile"], "legacy_nuclearid")
                self.assertEqual(maker["lambdaPairMassMode"], "fixed")
                self.assertEqual(float(maker["lambdaPairMass"]), 1.115683)
                self.assertEqual(maker["nChannels"], "4")
                mu = float(maker["lambdaSignalMean"])
                self.assertEqual(float(maker["lambdaSignalNSigma"]), 3.0)
                self.assertGreater(float(maker["lambdaSignalSigma"]), 0.0)
                width = float(maker["lambdaSignalSigma"]) * float(maker["lambdaSignalNSigma"])
                outer = float(maker["lambdaSidebandOuterFactor"]) * width
                bounds = [(1.05, 1.25), (mu-width, mu+width), (mu-outer, mu-width), (mu+width, mu+outer)]
                for i, (low, high) in enumerate(bounds):
                    for key in ("name", "partA", "partB", "enabled", "doMixing", "signalMin", "signalMax", "normQMin", "normQMax"):
                        self.assertIn(f"channel_{i}_{key}", maker)
                    self.assertAlmostEqual(float(maker[f"channel_{i}_signalMin"]), low)
                    self.assertAlmostEqual(float(maker[f"channel_{i}_signalMax"]), high)

    def test_entry_lifecycle_and_loader_order(self):
        for short, species in SPECIES.items():
            macro = (ROOT / f"analysis/anaFemtoLambda_{short}.C").read_text()
            self.assertIn('#include "StMaker/StFemtoMaker/StFemtoMaker.h"', macro)
            self.assertIn('new StFemtoMaker("femto"', macro)
            self.assertIn('SetLambdaProvider(createFemtoLambdaKfProvider())', macro)
            self.assertIn(f'run.Prepare("{species}")', macro)
            self.assertNotIn('StFemtoMakerLambdaNuclear', macro)
            runner = (ROOT / f"analysis/run_anaFemtoLambda_{short}.C").read_text()
            libs = ["libStarAnaConfig.so", "libStRefMultCorr.so", "libKFParticle.so",
                    "libStKfParticleCommon.so", "libStCommon.so", "libStFemtoMaker.so"]
            self.assertEqual([runner.index(x) for x in libs], sorted(runner.index(x) for x in libs))
            self.assertIn('gSystem->Exit', runner)

    def test_strict_event_count_and_per_file_covariance_guard(self):
        source = (ROOT / "analysis/FemtoLambdaRunSupport.h").read_text()
        self.assertIn('requested > available', source)
        self.assertIn('completed != requested', source)
        self.assertIn('GetTreeNumber()', source)
        self.assertIn('"TrackCovMatrix"', source)
        self.assertIn('maker->SetProcessingSucceeded(status == 0)', source)
        self.assertNotIn('Getenv(', source)

    def test_shell_syntax_and_explicit_arguments(self):
        scripts = [ROOT / "script/run_femtoLambda_common.sh", ROOT / "script/femtoLambda_environment.sh"]
        for short in SPECIES:
            for prefix in ("run", "singularity_run"):
                script = ROOT / f"script/{prefix}_anaFemtoLambda_{short}.sh"
                scripts.append(script)
                self.assertEqual(subprocess.run(["bash", str(script)], capture_output=True).returncode, 2)
        for script in scripts:
            subprocess.run(["bash", "-n", str(script)], check=True)
        env = (ROOT / "script/femtoLambda_environment.sh").read_text()
        self.assertNotIn('source ', env)
        self.assertNotIn('.current_mainconf', env)

    def test_qa_display_and_write_guards(self):
        source = (ROOT / "common/macro/checkHistAnaFemtoLambda.C").read_text()
        self.assertLess(source.index("cf->SetMaximum();"), source.index("cf->GetMaximum()"))
        self.assertIn("if (me->GetBinContent(bin) <= 0) continue;", source)
        self.assertIn("pdfInfo.fSize <= 0", source)
        self.assertIn("const TString outputDirectory = gSystem->DirName(outputPdf);", source)

    def test_nuclear_and_mixing_configuration(self):
        cfg = flat(ROOT / "config/cuts/nuclearid/nuclearid_femtoLambda_legacy.yaml")
        self.assertEqual(cfg["requireBestSpecies"], "true")
        self.assertEqual(cfg["nuclearMinNHitsDedx"], "15")
        self.assertEqual(cfg["nuclearMinPt"], "0.1")
        self.assertEqual(cfg["maxPOverQ"], "2.5")
        self.assertNotIn("MeanLambda", cfg)  # Windows have one authoritative maker definition.
        cfg = flat(ROOT / "config/cuts/mixing/mixing_femtoLambda_legacy.yaml")
        self.assertEqual(cfg["mixingMode"], "bufferAll")
        self.assertEqual(cfg["mixBothDirections"], "true")
        self.assertEqual(cfg["nCentralityBins"], "9")
        self.assertEqual(cfg["vzOutOfRangePolicy"], "clamp")
        self.assertLessEqual(int(cfg["maxMixEvents"]), int(cfg["bufferSize"]))


if __name__ == "__main__":
    unittest.main()
