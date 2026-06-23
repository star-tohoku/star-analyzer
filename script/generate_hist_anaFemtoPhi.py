#!/usr/bin/env python3
"""Generate config/hist/hist_anaFemtoPhi.yaml from existing per-species templates."""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
HIST_DIR = ROOT / "config" / "hist"

COMMON_STOP_KEYS = {
    "hHe4_Pt_PreFemtoCut",
    "hP_Pt_PreFemtoCut",
    "hDeuteron_Pt_PreFemtoCut",
}

BACHELOR_PREFIXES = ("hP_", "hDeuteron_", "hTriton_", "hHe3_", "hHe4_")
PHI_FEMTO_KEYS = {
    "hPhi_MKK_allForFemto",
    "hPhi_MKK_signal",
    "hPhi_MKK_leftSB",
    "hPhi_MKK_rightSB",
    "hPhi_MKK_rot",
    "hPhi_PtVsY_PostCut",
    "hPhi_PtVsY_signal",
    "hPhi_PtVsY_rot",
    "hPhiRot_MKK",
    "hPhiRot_Pt",
    "hPhiRot_Rapidity",
    "hPhiRot_NCand",
    "hPhiRot_DeltaPhiApplied",
}
SKIP_DUPLICATE = set(PHI_FEMTO_KEYS)


def parse_axes(text: str) -> str:
    m = re.search(r"^axes:\s*\n", text, re.MULTILINE)
    if not m:
        return ""
    start = m.end()
    m2 = re.search(r"^histograms:\s*$", text[start:], re.MULTILINE)
    end = start + m2.start() if m2 else len(text)
    return text[start:end]


def merge_axes(*axis_sections: str) -> str:
    blocks: dict[str, list[str]] = {}
    order: list[str] = []
    for section in axis_sections:
        current: list[str] = []
        current_name = ""
        for line in section.splitlines():
            if line.strip() == "axes:":
                continue
            if line.strip().startswith("#"):
                continue
            m = re.match(r"^  ([A-Za-z0-9_]+):", line)
            if m:
                name = m.group(1)
                if name.startswith("h"):
                    continue
                if current_name and current_name not in blocks:
                    blocks[current_name] = current
                    order.append(current_name)
                current_name = name
                current = [line]
            elif current_name:
                current.append(line)
        if current_name and current_name not in blocks:
            blocks[current_name] = current
            order.append(current_name)
    lines: list[str] = ["axes:"]
    for name in order:
        lines.extend(blocks[name])
    return "\n".join(lines) + "\n\n"


def parse_hist_blocks(text: str) -> tuple[str, dict[str, str]]:
    """Return (header through 'histograms:' line, {name: block})."""
    m = re.search(r"^histograms:\s*$", text, re.MULTILINE)
    if not m:
        raise RuntimeError("histograms: section not found")
    header = text[: m.end() + 1]
    body = text[m.end() :]
    blocks: dict[str, str] = {}
    for match in re.finditer(r"^  (h[A-Za-z0-9_]+):\s*\n", body, re.MULTILINE):
        name = match.group(1)
        start = match.start()
        nxt = re.search(r"^  h[A-Za-z0-9_]+:\s*\n", body[match.end() :], re.MULTILINE)
        end = match.end() + nxt.start() if nxt else len(body)
        blocks[name] = body[start:end].rstrip() + "\n"
    return header, blocks


def clone_he4_to_species(block: str, src_prefix: str, dst_prefix: str, replacements: list[tuple[str, str]]) -> str:
    out = block.replace(src_prefix, dst_prefix)
    for old, new in replacements:
        out = out.replace(old, new)
    return out


def he4_to_triton(block: str) -> str:
    return clone_he4_to_species(
        block,
        "hHe4_",
        "hTriton_",
        [
            ("NSigmaHe4", "NSigmaTriton"),
            ("n#sigma_{^{4}He}", "n#sigma_{^{3}H}"),
            ("^{4}He", "^{3}H"),
            ("N_{^{4}He}", "N_{^{3}H}"),
            ("hNHe4_vs_Cent9", "hNTriton_vs_Cent9"),
            ("Mass2He4", "Mass2He4"),  # reuse axis preset
        ],
    )


def he4_to_he3(block: str) -> str:
    return clone_he4_to_species(
        block,
        "hHe4_",
        "hHe3_",
        [
            ("NSigmaHe4", "NSigmaHe3"),
            ("n#sigma_{^{4}He}", "n#sigma_{^{3}He}"),
            ("^{4}He", "^{3}He"),
            ("N_{^{4}He}", "N_{^{3}He}"),
            ("hNHe4_vs_Cent9", "hNHe3_vs_Cent9"),
        ],
    )


def channel_block_from_template(tmpl: str, channel: str) -> str:
    return re.sub(r"phi_he4", channel, tmpl)


def extract_bachelor_blocks(blocks: dict[str, str], prefix: str) -> dict[str, str]:
    return {k: v for k, v in blocks.items() if k.startswith(prefix)}


SIGNAL_FEMTO_CHANNELS = [
    "phi_proton_signal",
    "phi_deuteron_signal",
    "phi_triton_signal",
    "phi_he3_signal",
    "phi_he4_signal",
]


def mkk_vs_kstar_block(channel: str, is_se: bool) -> str:
    kind = "SE" if is_se else "ME"
    name = f"hPhiMKK_vs_Kstar{kind}_{channel}"
    return f"""  {name}:
    xAxis: *MKK
    yAxis: *Kstar
    zAxis: *Cent9
    type: TH3F
    title: "#phi M_{{KK}} vs k^{{*}} {kind} vs cent9 ({channel});M_{{KK}} [GeV/c^{{2}}];k^{{*}} [GeV/c];cent9"
"""


def phi_pair_mom_angle_block(channel: str, tof_strict: bool) -> list[str]:
    suffix = "_tofStrict" if tof_strict else ""
    tag = "tofStrict" if tof_strict else "signal"
    h1d = f"hPhiPairMomAngle_{channel}{suffix}"
    h2d = f"hPhiPairMomAngle_vs_MKK_{channel}{suffix}"
    blocks = [
        f"""  {h1d}:
    axis: *OpeningAngle
    title: "#phi-bachelor momentum angle ({channel}, {tag});#theta_{{p}} [rad];Counts"
""",
        f"""  {h2d}:
    xAxis: *OpeningAngle
    yAxis: *MKK
    title: "#theta_{{p}} vs M_{{KK}} ({channel}, {tag});#theta_{{p}} [rad];M_{{KK}} [GeV/c^{{2}}]"
""",
    ]
    names = [h1d, h2d]
    return list(zip(names, blocks))


def main() -> None:
    he4_text = (HIST_DIR / "hist_anaFemtoPhi4He.yaml").read_text()
    proton_text = (HIST_DIR / "hist_anaFemtoPhiProton.yaml").read_text()
    deuteron_text = (HIST_DIR / "hist_anaFemtoPhiDeuteron.yaml").read_text()

    he4_header, he4_blocks = parse_hist_blocks(he4_text)
    _, proton_blocks = parse_hist_blocks(proton_text)
    _, deuteron_blocks = parse_hist_blocks(deuteron_text)

    merged_axes = merge_axes(
        parse_axes(he4_text),
        parse_axes(proton_text),
        parse_axes(deuteron_text),
    )

    common: dict[str, str] = {}
    for name, block in he4_blocks.items():
        if name.startswith("hHe4_") or name.startswith("hKstar") or name.startswith("hCF_"):
            continue
        if name in SKIP_DUPLICATE:
            continue
        if name.startswith("hNHe4_"):
            continue
        common[name] = block

    # Unified bachelor nσ vs p (plan §10.4)
    nsigma_extra = {
        "hNSigmaProtonVsP": """  hNSigmaProtonVsP:
    xAxis: *Momentum
    yAxis: *NSigmaY
    title: "n#sigma_{p} vs p (after ID cut);p [GeV/c];n#sigma_{p}"
""",
        "hNSigmaDeuteronVsP": """  hNSigmaDeuteronVsP:
    xAxis: *Momentum
    yAxis: *NSigmaY
    title: "n#sigma_{d} vs p (after ID cut);p [GeV/c];n#sigma_{d}"
""",
        "hNSigmaDeuteronVsP_All": """  hNSigmaDeuteronVsP_All:
    xAxis: *Momentum
    yAxis: *NSigmaY
    title: "n#sigma_{d} vs p (all tracks, no ID cut);p [GeV/c];n#sigma_{d}"
""",
        "hNSigmaTritonVsP": """  hNSigmaTritonVsP:
    xAxis: *Momentum
    yAxis: *NSigmaY
    title: "n#sigma_{^{3}H} vs p (after ID cut);p [GeV/c];n#sigma_{^{3}H}"
""",
        "hNSigmaTritonVsP_All": """  hNSigmaTritonVsP_All:
    xAxis: *Momentum
    yAxis: *NSigmaY
    title: "n#sigma_{^{3}H} vs p (all tracks, no ID cut);p [GeV/c];n#sigma_{^{3}H}"
""",
        "hNSigmaHe3VsP": """  hNSigmaHe3VsP:
    xAxis: *Momentum
    yAxis: *NSigmaY
    title: "n#sigma_{^{3}He} vs p (after ID cut);p [GeV/c];n#sigma_{^{3}He}"
""",
        "hNSigmaHe3VsP_All": """  hNSigmaHe3VsP_All:
    xAxis: *Momentum
    yAxis: *NSigmaY
    title: "n#sigma_{^{3}He} vs p (all tracks, no ID cut);p [GeV/c];n#sigma_{^{3}He}"
""",
    }
    for k, v in nsigma_extra.items():
        if k not in common:
            common[k] = v

    bachelor: dict[str, str] = {}
    bachelor.update(extract_bachelor_blocks(proton_blocks, "hP_"))
    bachelor["hNProton_vs_Cent9"] = proton_blocks.get(
        "hNProton_vs_Cent9",
        """  hNProton_vs_Cent9:
    xAxis: *Cent9
    yAxis: *NKaon
    title: "N_{p} vs cent9;cent9;N_{p}"
""",
    )
    bachelor.update(extract_bachelor_blocks(deuteron_blocks, "hDeuteron_"))
    bachelor["hNDeuteron_vs_Cent9"] = deuteron_blocks.get("hNDeuteron_vs_Cent9", he4_blocks.get("hNDeuteron_vs_Cent9", ""))
    bachelor.update(extract_bachelor_blocks(he4_blocks, "hHe4_"))
    bachelor["hNHe4_vs_Cent9"] = he4_blocks["hNHe4_vs_Cent9"]

    he4_bachelor = extract_bachelor_blocks(he4_blocks, "hHe4_")
    for name, block in he4_bachelor.items():
        tname = name.replace("hHe4_", "hTriton_")
        bachelor[tname] = he4_to_triton(block)
    bachelor["hNTriton_vs_Cent9"] = clone_he4_to_species(
        he4_blocks["hNHe4_vs_Cent9"], "hNHe4_", "hNTriton_", [("^{4}He", "^{3}H"), ("N_{^{4}He}", "N_{^{3}H}")]
    )
    for name, block in he4_bachelor.items():
        h3name = name.replace("hHe4_", "hHe3_")
        bachelor[h3name] = he4_to_he3(block)
    bachelor["hNHe3_vs_Cent9"] = clone_he4_to_species(
        he4_blocks["hNHe4_vs_Cent9"], "hNHe4_", "hNHe3_", [("^{4}He", "^{3}He"), ("N_{^{4}He}", "N_{^{3}He}")]
    )

    phi_femto: dict[str, str] = {}
    for key in sorted(PHI_FEMTO_KEYS):
        if key in he4_blocks:
            phi_femto[key] = he4_blocks[key]

    channels = [
        "phi_proton",
        "phi_proton_signal",
        "phi_proton_leftSB",
        "phi_proton_rightSB",
        "phi_rot_proton",
        "phi_deuteron",
        "phi_deuteron_signal",
        "phi_deuteron_leftSB",
        "phi_deuteron_rightSB",
        "phi_triton",
        "phi_triton_signal",
        "phi_triton_leftSB",
        "phi_triton_rightSB",
        "phi_he3",
        "phi_he3_signal",
        "phi_he3_leftSB",
        "phi_he3_rightSB",
        "phi_he4",
        "phi_he4_signal",
        "phi_he4_leftSB",
        "phi_he4_rightSB",
    ]

    channel_blocks: dict[str, str] = {}
    he4_channel_names = [k for k in he4_blocks if k.startswith("hKstar") or k.startswith("hCF_")]
    for ch in channels:
        for hk in he4_channel_names:
            suffix = hk.split("_phi_he4", 1)[-1] if "_phi_he4" in hk else ""
            base = hk.replace("phi_he4", "CHANNEL").replace("phi_rot_he4", "CHANNEL_ROT")
            if "CHANNEL_ROT" in base:
                if ch == "phi_rot_proton":
                    new_name = base.replace("CHANNEL_ROT", "phi_rot_proton")
                    if hk in he4_blocks:
                        channel_blocks[new_name] = channel_block_from_template(he4_blocks[hk], "phi_rot_proton")
                continue
            if "CHANNEL" not in base:
                continue
            new_name = base.replace("CHANNEL", ch)
            src_key = hk
            if hk not in he4_blocks:
                continue
            if ch.startswith("phi_proton"):
                src_key = hk.replace("phi_he4", "phi_proton")
                if src_key not in proton_blocks:
                    src_key = hk
                else:
                    channel_blocks[new_name] = channel_block_from_template(proton_blocks[src_key], ch)
                    continue
            if ch.startswith("phi_deuteron"):
                src_key = hk.replace("phi_he4", "phi_deuteron")
                if src_key in deuteron_blocks:
                    channel_blocks[new_name] = channel_block_from_template(deuteron_blocks[src_key], ch)
                    continue
            channel_blocks[new_name] = channel_block_from_template(he4_blocks[hk], ch)

    mkk_kstar_blocks: dict[str, str] = {}
    for ch in SIGNAL_FEMTO_CHANNELS:
        for is_se in (True, False):
            block = mkk_vs_kstar_block(ch, is_se)
            m = re.search(r"^  (h[A-Za-z0-9_]+):", block, re.MULTILINE)
            if m:
                mkk_kstar_blocks[m.group(1)] = block

    pair_angle_blocks: dict[str, str] = {}
    for ch in SIGNAL_FEMTO_CHANNELS:
        for tof_strict in (False, True):
            for name, block in phi_pair_mom_angle_block(ch, tof_strict):
                pair_angle_blocks[name] = block

    out_header = "# StFemtoMaker histogram definitions (auau3p85fxt_anaFemtoPhi unified)\n"
    out_header += "# 1D: axis: *Preset or nBins/min/max; title required\n"
    out_header += "# 2D: xAxis: *Preset, yAxis: *Preset; title required\n\n"
    out_header += merged_axes

    ordered_names: list[str] = []
    ordered_names.extend(sorted(common.keys()))
    ordered_names.extend(sorted(k for k in bachelor if k.startswith("hP_")))
    ordered_names.append("hNProton_vs_Cent9")
    for pref in ("hDeuteron_", "hTriton_", "hHe3_", "hHe4_"):
        ordered_names.extend(sorted(k for k in bachelor if k.startswith(pref)))
    for nc in ("hNDeuteron_vs_Cent9", "hNTriton_vs_Cent9", "hNHe3_vs_Cent9", "hNHe4_vs_Cent9"):
        if nc in bachelor and nc not in ordered_names:
            ordered_names.append(nc)
    ordered_names.extend(sorted(phi_femto.keys()))
    for ch in channels:
        for kind in ("hKstarSE_", "hKstarME_", "hCF_", "hKstarSEVsCent_", "hKstarMEVsCent_"):
            nm = kind + ch
            if nm in channel_blocks:
                ordered_names.append(nm)
    for ch in SIGNAL_FEMTO_CHANNELS:
        for kind in ("hPhiMKK_vs_KstarSE_", "hPhiMKK_vs_KstarME_"):
            nm = kind + ch
            if nm in mkk_kstar_blocks:
                ordered_names.append(nm)
    for ch in SIGNAL_FEMTO_CHANNELS:
        for tof_strict in (False, True):
            for name, _ in phi_pair_mom_angle_block(ch, tof_strict):
                if name in pair_angle_blocks:
                    ordered_names.append(name)

    merged: dict[str, str] = {}
    merged.update(common)
    merged.update(bachelor)
    merged.update(phi_femto)
    merged.update(channel_blocks)
    merged.update(mkk_kstar_blocks)
    merged.update(pair_angle_blocks)

    lines = [out_header.rstrip(), "", "histograms:"]
    for name in ordered_names:
        if name in merged:
            lines.append(merged[name].rstrip())
            lines.append("")
    out_path = HIST_DIR / "hist_anaFemtoPhi.yaml"
    out_text = "\n".join(lines).rstrip() + "\n"
    if "\nhistograms:\n" not in out_text and not out_text.startswith("histograms:\n"):
        raise RuntimeError("generated YAML missing top-level histograms: section")
    out_path.write_text(out_text)
    print(f"Wrote {out_path} ({len(ordered_names)} histograms)")


if __name__ == "__main__":
    main()
