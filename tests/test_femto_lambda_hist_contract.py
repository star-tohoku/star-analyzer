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
        top = yaml.safe_load((ROOT / ("config/hist/hist_femtoLambda_" + species + "_kf.yaml")).read_text())["histograms"]
        nuclear = yaml.safe_load((ROOT / ("config/hist/hist_femtoLambda_" + species + "_nuclear.yaml")).read_text())["histograms"]
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
        assert set(nuclear) == legacy_expected | {"hDphiDeta_proton_" + legacy + "_EventField"}
        for mode in ("SE", "ME"):
            for region in ("", "_signal", "_leftSB", "_rightSB"):
                channel = "lambda_" + species + region
                assert "hKstar" + mode + "_" + channel in top
                assert "hKstar" + mode + "VsCent_" + channel in top
        assert len(top) == 67
        assert len(nuclear) == 113
    print("PASS original Lambda histogram contract:", tested, "legacy key/axis/class comparisons across four species")


if __name__ == "__main__":
    main()

