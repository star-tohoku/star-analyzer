#!/usr/bin/env python3
"""Mechanically import the coherent STAR KFParticle reconstruction core.

External headers stay in the global namespace. Each upstream C++ declaration
is enclosed in star_analyzer_kfp, including bundled SIMD and allocator types.
Default physics selections and STAR conditional branches are preserved. Two
explicit Lambda-only study controls are added reproducibly; see the generated
manifest and StRoot/KFParticle/PROVENANCE.md for snapshot details.
"""

import argparse
import datetime
import hashlib
import json
from pathlib import Path
import re
import sys


CORE = (
    "KFParticle", "KFPTrack", "KFPVertex", "KFParticleDatabase", "KFVertex",
    "KFPTrackVector", "KFPEmcCluster", "KFParticleSIMD",
    "KFParticlePVReconstructor", "KFParticleFinder", "KFParticleTopoReconstructor",
)
HEADERS = ("KFParticleDef.h", "KFParticleMath.h", "KFParticleField.h",
           "KFPInputData.h", "KFPSimdAllocator.h")
NAMESPACE = "star_analyzer_kfp"
OPEN = "namespace " + NAMESPACE + " { // local ABI isolation\n"
CLOSE = "} // namespace " + NAMESPACE + "\n"


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def lambda_study_controls(name, text):
    """Apply only audited Lambda selection hooks, failing on upstream drift."""
    def replace(old, new):
        nonlocal text
        if text.count(old) != 1:
            raise RuntimeError("Lambda study transform anchor changed: " + name + ": " + old)
        text = text.replace(old, new)

    if name == "KFParticleTopoReconstructor.h":
        replace("  ~KFParticleTopoReconstructor();",
                "  ~KFParticleTopoReconstructor();\n"
                "  void SetLambdaTopoChi2NdfCut(float cut) { fLambdaTopoChi2NdfCut = cut; }\n"
                "  float GetLambdaTopoChi2NdfCut() const { return fLambdaTopoChi2NdfCut; }")
        replace("    fNThreads = a.fNThreads;",
                "    fNThreads = a.fNThreads;\n    fLambdaTopoChi2NdfCut = a.fLambdaTopoChi2NdfCut;")
        replace("  void CopyCuts(const KFParticleTopoReconstructor* topo) { fKFParticleFinder->CopyCuts(topo->fKFParticleFinder); }",
                "  void CopyCuts(const KFParticleTopoReconstructor* topo) { fKFParticleFinder->CopyCuts(topo->fKFParticleFinder); fLambdaTopoChi2NdfCut = topo->fLambdaTopoChi2NdfCut; }")
        replace("}__attribute__((aligned(sizeof(float32_v)))); // class KFParticleTopoReconstructor",
                "  // Appended local study setting; the upstream default is unchanged.\n"
                "  float fLambdaTopoChi2NdfCut = 3.f;\n"
                "}__attribute__((aligned(sizeof(float32_v)))); // class KFParticleTopoReconstructor")
        replace("  {\n  }\n  \n  /** Copy cuts from KF Particle Finder",
                "  {\n    fLambdaTopoChi2NdfCut = a.fLambdaTopoChi2NdfCut;\n  }\n  \n  /** Copy cuts from KF Particle Finder")
    elif name == "KFParticleTopoReconstructor.cxx":
        replace("      if(tmp.Chi2()/tmp.NDF()<3.)",
                "      const float topoCut = abs(fParticles[iParticle].GetPDG()) == 3122\n"
                "          ? fLambdaTopoChi2NdfCut : 3.f;\n"
                "      if(tmp.Chi2()/tmp.NDF()<topoCut)")
    elif name == "KFParticleFinder.h":
        replace("  void SetLCut(float cut) { fLCut = cut; }",
                "  void SetApplyLambdaGeometryCuts(bool apply) { fApplyLambdaGeometryCuts = apply; }\n"
                "  bool GetApplyLambdaGeometryCuts() const { return fApplyLambdaGeometryCuts; }\n"
                "  void SetLCut(float cut) { fLCut = cut; }")
        replace("    fLCut = finder->fLCut;",
                "    fLCut = finder->fLCut;\n    fApplyLambdaGeometryCuts = finder->fApplyLambdaGeometryCuts;")
        replace("  KFParticleFinder(const KFParticleFinder&);",
                "  // Appended local study switch: only Lambda/anti-Lambda geometry gates.\n"
                "  bool fApplyLambdaGeometryCuts = true;\n\n"
                "  KFParticleFinder(const KFParticleFinder&);")
    elif name == "KFParticleFinder.cxx":
        # Restrict these three gates to the two-daughter ConstructV0 path.
        begin = text.index("inline void KFParticleFinder::ConstructV0(")
        end = text.index("inline void KFParticleFinder::SaveV0PrimSecCand(", begin)
        prefix, suffix = text[:begin], text[end:]
        text = text[begin:end]
        replace("  saveParticle &= (lMin < 200.f);",
                "  const mask32_v skipLambdaGeometry = fApplyLambdaGeometryCuts\n"
                "      ? (int32_v(0) == int32_v(1)) : (abs(mother.PDG()) == int32_v(3122));\n"
                "  saveParticle &= (lMin < 200.f) || skipLambdaGeometry;")
        replace("  saveParticle &= ((!isPrimary) && isParticleFromVertex) || isPrimary;",
                "  saveParticle &= ((!isPrimary) && isParticleFromVertex) || isPrimary || skipLambdaGeometry;")
        replace("  saveParticle &= ( ((isK0 || isLambda || isHyperNuclei) && lMin > float32_v(fLCut)) || !(isK0 || isLambda || isHyperNuclei) );",
                "  saveParticle &= ( ((isK0 || isLambda || isHyperNuclei) && lMin > float32_v(fLCut)) || !(isK0 || isLambda || isHyperNuclei) ) || skipLambdaGeometry;")
        text = prefix + text + suffix
        replace("                  active[iPDGPos] &= (dr < float32_v(fDistanceCut));",
                "                  const mask32_v skipLambdaGeometry = fApplyLambdaGeometryCuts\n"
                "                      ? (int32_v(0) == int32_v(1)) : (abs(motherPDG) == int32_v(3122));\n"
                "                  active[iPDGPos] &= (dr < float32_v(fDistanceCut)) || skipLambdaGeometry;")
        replace("                  active[iPDGPos] &= (p1p2 > -p12);\n                  active[iPDGPos] &= (p1p2 > -p22);",
                "                  active[iPDGPos] &= (p1p2 > -p12) || skipLambdaGeometry;\n"
                "                  active[iPDGPos] &= (p1p2 > -p22) || skipLambdaGeometry;")
    return text


def isolate(name, data, guard_names):
    text = lambda_study_controls(name, data.decode("utf-8"))
    # Unlike the other internal KF types, upstream KFPTrack ignores the
    # standalone switch for ClassDef. Preserve TObject/STAR layout, but apply
    # the same dictionary opt-out as KFParticle/KFVertex.
    if name == "KFPTrack.h":
        old = "#ifdef __ROOT__\n  ClassDef(KFPTrack,1)\n#endif"
        new = "#if defined(__ROOT__) && !defined(KFParticleStandalone)\n  ClassDef(KFPTrack,1)\n#endif"
        if text.count(old) != 1:
            raise RuntimeError("KFPTrack dictionary guard changed upstream")
        text = text.replace(old, new)
    # Guard names must not collide with StarRoot's old global KFParticle.h.
    for name in guard_names:
        text = re.sub(r"\b" + re.escape(name) + r"\b",
                      "STAR_ANALYZER_VENDOR_" + name, text)
    lines = []
    for line in text.splitlines(keepends=True):
        # All includes in this audited snapshot occur at C++ global scope.
        # Close/reopen inside the SAME preprocessor branch, including disabled
        # optional integrations, so conditional compilation remains unchanged.
        if re.match(r"\s*#\s*include\b", line):
            lines.extend((CLOSE, line.rstrip("\n") + "\n", OPEN))
        else:
            lines.append(line)
    return ("// Generated by tools/kfparticle/import_core.py; see PROVENANCE.md.\n"
            + OPEN + "".join(lines).rstrip("\n") + "\n" + CLOSE).encode("utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path,
                        default=Path("/star/nfs4/AFS/star/packages/.DEV2/StRoot/KFParticle"))
    parser.add_argument("--snapshot-date", default="2026-09-07")
    parser.add_argument("--check", action="store_true",
                        help="verify generated files against source without writing")
    args = parser.parse_args()
    datetime.date.fromisoformat(args.snapshot_date)
    root = Path(__file__).resolve().parents[2]
    destination = root / "StRoot" / "KFParticle"
    source = args.source.resolve(strict=True)
    names = [name + ext for name in CORE for ext in (".h", ".cxx")]
    names.extend(HEADERS)
    names.extend(str(path.relative_to(source)) for path in
                 sorted((source / "KFPSimd").rglob("*.h")))
    names.extend(("KFPSimd/COPYING", "KFPSimd/compile_flags.txt"))
    contents = {name: (source / name).read_bytes() for name in names}
    guards = []
    for name, content in contents.items():
        if name.endswith(".h"):
            match = re.search(rb"^\s*#\s*ifndef\s+(\w+)", content, re.M)
            if not match:
                raise RuntimeError("Missing include guard: " + name)
            guards.append(match.group(1).decode("ascii"))
    generated = {name: (isolate(name, content, guards)
                        if name.endswith((".h", ".cxx")) else content)
                 for name, content in contents.items()}
    source_id = sha256("".join(name + ":" + sha256(contents[name]) + "\n"
                               for name in sorted(contents)).encode("utf-8"))
    manifest_path = destination / "SNAPSHOT.json"
    previous = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    previous_files = previous.get("files", {})
    for name, data in generated.items():
        target = destination / name
        if not target.exists():
            if args.check:
                raise RuntimeError("Missing vendor file: " + name)
            continue
        current = target.read_bytes()
        if args.check:
            if current != data:
                raise RuntimeError("Vendor file differs from mechanical import: " + name)
        elif current != data and current != contents[name]:
            expected = previous_files.get(name, {}).get("vendored_sha256")
            if expected != sha256(current):
                raise RuntimeError("Refusing to overwrite unrecorded local edit: " + name)
    manifest = {
        "source": str(source), "snapshot_date": args.snapshot_date,
        "source_id": source_id, "namespace": NAMESPACE,
        "transform": "namespace-guards-kfptrack-dictionary-lambda-study-controls-v3",
        "build_defines": ["__ROOT__", "KFParticleStandalone", "HomogeneousField"],
        "simd": "bundled KFPSimd SSE4.1; -msse4.1, not -march=native",
        "translation_units": [name + ".cxx" for name in CORE],
        "files": {name: {"upstream_sha256": sha256(contents[name]),
                         "vendored_sha256": sha256(data)}
                  for name, data in sorted(generated.items())},
    }
    manifest_bytes = (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode("utf-8")
    fingerprint = ("// Generated by tools/kfparticle/import_core.py\n"
                   "#ifndef STAR_ANALYZER_KFP_SNAPSHOT_H\n"
                   "#define STAR_ANALYZER_KFP_SNAPSHOT_H\n"
                   '#define STAR_ANALYZER_KFP_SOURCE_ID "' + source_id + '"\n'
                   '#define STAR_ANALYZER_KFP_SNAPSHOT_DATE "' + args.snapshot_date + '"\n'
                   "#endif\n").encode("utf-8")
    fingerprint_path = destination / "KFParticleSnapshot.h"
    if args.check:
        if not manifest_path.exists() or manifest_path.read_bytes() != manifest_bytes:
            raise RuntimeError("Snapshot manifest does not match source and transform")
        if not fingerprint_path.exists() or fingerprint_path.read_bytes() != fingerprint:
            raise RuntimeError("Snapshot fingerprint header mismatch")
        print("Verified %d files, %d translation units, namespace %s" %
              (len(generated), len(CORE), NAMESPACE))
        return
    for name, data in generated.items():
        target = destination / name
        target.parent.mkdir(parents=True, exist_ok=True)
        if not target.exists() or target.read_bytes() != data:
            target.write_bytes(data)
    manifest_path.write_bytes(manifest_bytes)
    fingerprint_path.write_bytes(fingerprint)
    print("Imported %d files, %d translation units, namespace %s" %
          (len(generated), len(CORE), NAMESPACE))


if __name__ == "__main__":
    try:
        main()
    except (OSError, RuntimeError, ValueError) as error:
        sys.exit("KFParticle import: " + str(error))
