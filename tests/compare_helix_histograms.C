// Read-only ROOT 5/CINT comparator. Run on equivalent Helix control samples:
// root4star -l -b -q 'tests/compare_helix_histograms.C("before.root","after.root")'
// Exact comparison, not a statistical equivalence test. Differences caused by
// independently randomized input/corrections are still reported as differences.
#include "TFile.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TClass.h"
#include "TH1.h"
#include "TAxis.h"
#include "TList.h"
#include "TObjString.h"
#include "TMath.h"
#include "TSystem.h"
#include "TString.h"
#include <iostream>
#include <cstring>

void HelixCompareDifference(const TString& path, const char* quantity,
                            Double_t before, Double_t after, Long64_t& differences) {
  ++differences;
  // Bound noisy output, but keep checking every bin and count every difference.
  if (differences <= 30)
    std::cerr << "DIFF " << path << " " << quantity << ": " << before << " vs " << after << std::endl;
}

Bool_t HelixCompareNumber(Double_t before, Double_t after) {
  return TMath::Finite(before) && TMath::Finite(after) && before == after;
}

void HelixCompareAxis(const TAxis* before, const TAxis* after,
                       const TString& path, const char* name, Long64_t& differences) {
  if (before->GetNbins() != after->GetNbins()) {
    HelixCompareDifference(path, Form("%s bins", name), before->GetNbins(), after->GetNbins(), differences);
    return;
  }
  const Int_t bins = before->GetNbins();
  for (Int_t bin = 1; bin <= bins + 1; ++bin) {
    const Double_t a = before->GetBinLowEdge(bin);
    const Double_t b = after->GetBinLowEdge(bin);
    if (!HelixCompareNumber(a, b))
      HelixCompareDifference(path, Form("%s edge %d", name, bin), a, b, differences);
    if (bin <= bins && std::strcmp(before->GetBinLabel(bin), after->GetBinLabel(bin)) != 0) {
      ++differences;
      if (differences <= 30)
        std::cerr << "DIFF " << path << " " << name << " label " << bin << std::endl;
    }
  }
}

// ROOT 5.34 has no public GetNcells(). Its public GetBin maps the final
// overflow indices to the highest global bin of TH1/TH2/TH3 histograms.
Int_t HelixCompareCells(const TH1* histogram) {
  return histogram->GetBin(histogram->GetNbinsX() + 1,
      histogram->GetDimension() >= 2 ? histogram->GetNbinsY() + 1 : 0,
      histogram->GetDimension() >= 3 ? histogram->GetNbinsZ() + 1 : 0) + 1;
}

void HelixCompareHistogram(const TH1* before, const TH1* after,
                            const TString& path, Long64_t& differences, Long64_t& cells) {
  if (std::strcmp(before->ClassName(), after->ClassName()) != 0) {
    ++differences;
    if (differences <= 30)
      std::cerr << "DIFF " << path << " class: " << before->ClassName()
                << " vs " << after->ClassName() << std::endl;
  }
  if (before->GetDimension() != after->GetDimension()) {
    HelixCompareDifference(path, "dimensions", before->GetDimension(), after->GetDimension(), differences);
    return;
  }
  HelixCompareAxis(before->GetXaxis(), after->GetXaxis(), path, "x", differences);
  if (before->GetDimension() >= 2)
    HelixCompareAxis(before->GetYaxis(), after->GetYaxis(), path, "y", differences);
  if (before->GetDimension() >= 3)
    HelixCompareAxis(before->GetZaxis(), after->GetZaxis(), path, "z", differences);
  if (!HelixCompareNumber(before->GetEntries(), after->GetEntries()))
    HelixCompareDifference(path, "entries", before->GetEntries(), after->GetEntries(), differences);
  if (HelixCompareCells(before) != HelixCompareCells(after)) {
    HelixCompareDifference(path, "global cells", HelixCompareCells(before), HelixCompareCells(after), differences);
    return;
  }
  // The global-cell range covers every dimension's underflow and overflow,
  // including mixed under/overflow corners of TH2/TH3.
  for (Int_t cell = 0; cell < HelixCompareCells(before); ++cell) {
    ++cells;
    const Double_t a = before->GetBinContent(cell);
    const Double_t b = after->GetBinContent(cell);
    const Double_t ea = before->GetBinError(cell);
    const Double_t eb = after->GetBinError(cell);
    if (!HelixCompareNumber(a, b))
      HelixCompareDifference(path, Form("cell %d content", cell), a, b, differences);
    if (!HelixCompareNumber(ea, eb))
      HelixCompareDifference(path, Form("cell %d error", cell), ea, eb, differences);
  }
}

void HelixCompareDirectory(TDirectory* before, TDirectory* after, const TString& prefix,
                            Long64_t& histograms, Long64_t& nonhistBefore,
                            Long64_t& nonhistAfter, Long64_t& differences, Long64_t& cells) {
  TList seenBefore;
  seenBefore.SetOwner(kTRUE);
  TIter nextBefore(before->GetListOfKeys());
  TKey* key = 0;
  while ((key = (TKey*)nextBefore())) {
    if (seenBefore.FindObject(key->GetName())) continue; // Compare latest key cycle once.
    seenBefore.Add(new TObjString(key->GetName()));
    const TString path = prefix + key->GetName();
    TObject* a = before->Get(key->GetName());
    if (!a) { ++differences; std::cerr << "DIFF unreadable reference key " << path << std::endl; continue; }
    if (a->InheritsFrom(TDirectory::Class())) {
      TObject* b = after->Get(key->GetName());
      if (!b || !b->InheritsFrom(TDirectory::Class())) {
        ++differences; std::cerr << "DIFF missing candidate directory " << path << std::endl;
      } else {
        HelixCompareDirectory((TDirectory*)a, (TDirectory*)b, path + "/", histograms,
                                nonhistBefore, nonhistAfter, differences, cells);
      }
    } else if (a->InheritsFrom(TH1::Class())) {
      ++histograms;
      TObject* b = after->Get(key->GetName());
      if (!b || !b->InheritsFrom(TH1::Class())) {
        ++differences; std::cerr << "DIFF missing candidate histogram " << path << std::endl;
      } else {
        HelixCompareHistogram((TH1*)a, (TH1*)b, path, differences, cells);
      }
    } else {
      ++nonhistBefore; // Metadata/trees are counted but deliberately not compared.
    }
  }
  TList seenAfter;
  seenAfter.SetOwner(kTRUE);
  TIter nextAfter(after->GetListOfKeys());
  while ((key = (TKey*)nextAfter())) {
    if (seenAfter.FindObject(key->GetName())) continue;
    seenAfter.Add(new TObjString(key->GetName()));
    TObject* b = after->Get(key->GetName());
    if (!b) { ++differences; std::cerr << "DIFF unreadable candidate key " << prefix << key->GetName() << std::endl; continue; }
    if (b->InheritsFrom(TH1::Class())) {
      TObject* a = before->Get(key->GetName());
      if (!a || !a->InheritsFrom(TH1::Class())) {
        ++differences; std::cerr << "DIFF extra candidate histogram " << prefix << key->GetName() << std::endl;
      }
    } else if (b->InheritsFrom(TDirectory::Class())) {
      TObject* a = before->Get(key->GetName());
      if (!a || !a->InheritsFrom(TDirectory::Class())) {
        ++differences; std::cerr << "DIFF extra candidate directory " << prefix << key->GetName() << std::endl;
      }
    } else {
      ++nonhistAfter;
    }
  }
}

void compare_helix_histograms(const char* beforePath, const char* afterPath) {
  TFile* before = TFile::Open(beforePath, "READ");
  TFile* after = TFile::Open(afterPath, "READ");
  if (!before || before->IsZombie() || !after || after->IsZombie()) {
    std::cerr << "FAIL Helix histogram comparison: cannot open both input files" << std::endl;
    delete before; delete after;
    gSystem->Exit(2); return;
  }
  Long64_t histograms = 0, nonhistBefore = 0, nonhistAfter = 0, differences = 0, cells = 0;
  HelixCompareDirectory(before, after, "", histograms, nonhistBefore, nonhistAfter, differences, cells);
  if (!histograms) { ++differences; std::cerr << "DIFF no reference histograms" << std::endl; }
  std::cout << "Helix histogram comparison: reference=" << beforePath << " candidate=" << afterPath << std::endl;
  std::cout << "histograms=" << histograms << " cells=" << cells << " differences=" << differences
            << " nonhistKeysNotCompared=" << nonhistBefore << "/" << nonhistAfter << std::endl;
  std::cout << (differences ? "FAIL" : "PASS")
            << " exact TH1-subclass classes/dimensions/axes/labels/entries/contents/errors, including under/overflow"
            << std::endl;
  before->Close(); after->Close(); delete before; delete after;
  gSystem->Exit(differences ? 1 : 0);
}

