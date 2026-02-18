---
name: stmaker-add-histograms
description: Steps to add new histograms to StPhiMaker (or similar StMaker) in this STAR analysis repo. Use when adding histograms, new axes, or checkHist pages.
---

# Adding Histograms (StPhiMaker / StMaker)

## Steps

1. **config/hist/hist_auau19_anaPhi.yaml**
   - Add axis under `axes:` if needed (e.g. `NKaon` with bins).
   - Add entry under `histograms:` with `axis:` (1D) or `xAxis`/`yAxis` (2D), and `title`.

2. **StMaker/StPhiMaker/StPhiMaker.cxx**
   - Call `m_histManager->Fill("HistName", ...)` in the right place:
     - Event-level (e.g. N(K+), N(K-)): after building the kaon list, before/after pair loop.
     - Pair-level (e.g. M_KK in rapidity/mult bins): inside the KK pair loop; use existing patterns for rapidity windows or multiplicity bins.

3. **common/macro/checkHistAnaPhi.C**
   - Add a page or panels: `fin->Get("HistName")` and draw. Use `Draw("colz")` for 2D; set `gPad->SetLogy(0)` when log scale is not wanted.

## Conventions

- Histogram key in YAML and in `Fill("...")` / `Get("...")` must match.
- Rebuild after changing YAML or .cxx; macro changes need no build.
