#!/usr/bin/env python3
"""Verify immutable legacy code/configuration against the pre-implementation archive."""
import hashlib
import json
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT.parent / "implementation_backups/femto_lambda_20260930/before.tar.gz"
PREFIXES = (
    "config/", "StMaker/StLambdaMaker/", "StMaker/StNuclearIdMaker/",
    "StMaker/StLambdaNuclearMixMaker/", "StMaker/StLambdaNuclearMaker/",
)
EXACT = {
    "analysis/anaLambdaNuclearId.C", "analysis/run_anaLambdaNuclearId.C",
    "analysis/anaLambda.C", "analysis/run_anaLambda.C",
    "StMaker/common/StLambdaV0Reconstruction.h",
    "StMaker/common/StLambdaV0Reconstruction.cxx",
    "StMaker/kfparticle/StPicoKFParticleInterface.cxx",
    "include/cuts/KfParticleCutConfig.h", "src/cuts/KfParticleCutConfig.cpp",
}
checked = []
with tarfile.open(ARCHIVE, "r:gz") as archive:
    for member in archive:
        name = member.name
        if not member.isfile() or not (name.startswith(PREFIXES) or name in EXACT):
            continue
        before = archive.extractfile(member).read()
        after = (ROOT / name).read_bytes()
        if before != after:
            raise SystemExit("FAIL: protected file changed: " + name)
        checked.append({"path": name, "sha256": hashlib.sha256(after).hexdigest()})
print(json.dumps({"status": "PASS", "protected_file_count": len(checked),
    "archive_sha256": hashlib.sha256(ARCHIVE.read_bytes()).hexdigest(), "files": checked}, indent=2))
