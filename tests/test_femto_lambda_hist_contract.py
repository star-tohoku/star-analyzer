#!/usr/bin/env python3
"""Read-only contract check against ORIGINAL anaLambdaNuclearId histogram YAMLs."""
import json
import pathlib
import re
import yaml

ROOT = pathlib.Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "mdfiles/femto/plans/femto_lambda_legacy_histogram_manifest_20260930.json"


def signature(spec):
    axes = [spec[k] for k in ("axis", "xAxis", "yAxis", "zAxis") if k in spec]
    if not axes:
        axes = [{k: spec[k] for k in ("nBins", "min", "max")}]
    return spec.get("type", "TH2F" if len(axes) == 2 else "TH1F"), axes


def main():
    manifest = json.loads(MANIFEST.read_text())
    assert "anaLambdaNuclearId.C" in manifest["reference"]
    old = manifest["histograms"]
    for source in {item["source"] for item in old}:
        configured = yaml.safe_load((ROOT / source).read_text())["histograms"]
        assert set(configured) == {item["key"] for item in old if item["source"] == source}
        for item in (v for v in old if v["source"] == source):
            assert signature(configured[item["key"]]) == (item["type"], item["axes"]), item["key"]

    tested = 0
    for species, legacy in (("deuteron", "d"), ("triton", "t"), ("he3", "3He"), ("he4", "4He")):
        top_doc = yaml.safe_load((ROOT / ("config/hist/hist_femtoLambda_" + species + "_kf.yaml")).read_text())
        top = top_doc["histograms"]
        nuclear_doc = yaml.safe_load((ROOT / ("config/hist/hist_femtoLambda_" + species + "_nuclear.yaml")).read_text())
        nuclear = nuclear_doc["histograms"]
        legacy_expected = set()
        for item in old:
            key = item["key"]
            match = re.match(r"^(?:hKstar(?:Mass)?|hQlab)_(?:Mixed_)?(d|t|3He|4He)(?:_|$)", key)
            if match is None:
                match = re.match(r"^hDphiDeta_proton_(d|t|3He|4He)$", key)
            if match and match.group(1) != legacy:
                continue
            dest = top if item["directory"] == "/" else nuclear
            assert key in dest, (species, key)
            assert signature(dest[key]) == (item["type"], item["axes"]), (species, key)
            if item["directory"] != "/":
                legacy_expected.add(key)
            tested += 1
        mt_expected = set()
        assert nuclear_doc["axes"]["PairMt"] == {"nBins": 200, "min": 0.0, "max": 10.0}
        for mixed in ("", "Mixed_"):
            for region in ("", "_SBPos", "_SBNeg", "_signalLow", "_signalHigh"):
                for cent in range(9):
                    key = "hKstarMt_" + mixed + legacy + region + "_CentBin" + str(cent)
                    mt_expected.add(key)
                    # ProjectionX must have exactly the original per-centrality k* axis.
                    original_region = "" if region in ("_signalLow", "_signalHigh") else region
                    original = "hKstar_" + mixed + legacy + original_region + "_CentBin" + str(cent)
                    assert signature(nuclear[key]) == ("TH2D", [nuclear[original]["axis"],
                                                               nuclear_doc["axes"]["PairMt"]]), (species, key)
        assert len(mt_expected) == 90
        assert set(nuclear) == legacy_expected | {"hDphiDeta_proton_" + legacy + "_EventField"} | mt_expected
        for mode in ("SE", "ME"):
            for region in ("", "_signal", "_leftSB", "_rightSB"):
                channel = "lambda_" + species + region
                assert "hKstar" + mode + "_" + channel in top
                assert "hKstar" + mode + "VsCent_" + channel in top
            for region in ("_signalLow", "_signalHigh"):
                for centrality in ("", "VsCent"):
                    stem = "hKstar" + mode + centrality + "_lambda_" + species
                    assert signature(top[stem + region]) == signature(top[stem + "_signal"]), (species, stem, region)
        acceptance = ("hLambda_PtVsYLab_signal", "hLambdaProton_PtVsYLab_signal",
                      "hLambdaPion_PtVsYLab_signal", "hNucleus_PtVsYLab_" + species)
        for key in acceptance:
            assert signature(top[key]) == ("TH2D", [top_doc["axes"]["YLab"], top_doc["axes"]["Pt"]]), key
        assert top_doc["axes"]["YLab"]["min"] < 0 < top_doc["axes"]["YLab"]["max"]
        assert len(top) == 79
        assert len(nuclear) == 203
        assert len(top) + len(nuclear) == 282
    print("PASS original Lambda histogram contract:", tested, "legacy key/axis/class comparisons across four species")
    print("PASS k*-mT histogram contract: 90 new TH2D keys per species, 282 total histograms per output")


if __name__ == "__main__":
    main()
