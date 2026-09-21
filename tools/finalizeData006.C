// Finalize merged DATA-006 downstream histograms into fitter-facing ROOT/CSV/QA products.
// ROOT 5 / C++98 compatible.
#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TMath.h"
#include "TNamed.h"
#include "TParameter.h"
#include "TString.h"
#include "TSystem.h"
#include "Data006Config.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

namespace {

Double_t QuietNaN() { return std::numeric_limits<Double_t>::quiet_NaN(); }

TH1D* ProjectCent(TH2D* h, const Data006CentralityConfig& c, const Char_t* name) {
  if (!h) return 0;
  return h->ProjectionX(name, c.cent9Min + 1, c.cent9Max + 1, "e");
}

TH1D* ProjectMt(TH3F* h, const Data006CentralityConfig& c, Double_t kmax,
                const Char_t* name) {
  if (!h) return 0;
  const Int_t ky1 = h->GetYaxis()->FindFixBin(0.0 + 1e-12);
  const Int_t ky2 = h->GetYaxis()->FindFixBin(kmax - 1e-12);
  return h->ProjectionX(name, ky1, ky2, c.cent9Min + 1, c.cent9Max + 1, "e");
}

TH2D* ProjectKstarMt(TH3F* h, const Data006CentralityConfig& c, Double_t kmax,
                     const Char_t* name) {
  if (!h) return 0;
  const Int_t nk = h->GetYaxis()->FindFixBin(kmax - 1e-12);
  TH2D* out = new TH2D(name, "k* vs pair m_{T}/2;k* (GeV/c);m_{T,pair}/2 (GeV/c^{2})",
                       nk, h->GetYaxis()->GetXmin(), h->GetYaxis()->GetBinUpEdge(nk),
                       h->GetXaxis()->GetNbins(), h->GetXaxis()->GetXmin(),
                       h->GetXaxis()->GetXmax());
  out->Sumw2();
  // The delivered 2D object is explicitly truncated at kmax. Fill its
  // underflow and in-range bins only; do not mislabel source bin nk+1 as the
  // aggregate overflow of the truncated axis.
  for (Int_t ik = 0; ik <= nk; ++ik) {
    for (Int_t im = 0; im <= h->GetXaxis()->GetNbins() + 1; ++im) {
      Double_t sum = 0.0, err2 = 0.0;
      for (Int_t iz = c.cent9Min + 1; iz <= c.cent9Max + 1; ++iz) {
        const Double_t x = h->GetBinContent(im, ik, iz);
        const Double_t e = h->GetBinError(im, ik, iz);
        sum += x;
        err2 += e * e;
      }
      out->SetBinContent(ik, im, sum);
      out->SetBinError(ik, im, TMath::Sqrt(err2));
    }
  }
  return out;
}

TH1D* ProjectPairY(TH2D* h, const Data006CentralityConfig& c, const Char_t* name) {
  if (!h) return 0;
  return h->ProjectionX(name, c.cent9Min + 1, c.cent9Max + 1, "e");
}

TH2D* ProjectXYCent(TH3F* h, const Data006CentralityConfig& c, const Char_t* name,
                    const Char_t* title) {
  if (!h) return 0;
  TH2D* out = new TH2D(name, title, h->GetXaxis()->GetNbins(), h->GetXaxis()->GetXmin(),
                       h->GetXaxis()->GetXmax(), h->GetYaxis()->GetNbins(),
                       h->GetYaxis()->GetXmin(), h->GetYaxis()->GetXmax());
  out->Sumw2();
  for (Int_t ix = 0; ix <= h->GetXaxis()->GetNbins() + 1; ++ix) {
    for (Int_t iy = 0; iy <= h->GetYaxis()->GetNbins() + 1; ++iy) {
      Double_t sum = 0.0, err2 = 0.0;
      for (Int_t iz = c.cent9Min + 1; iz <= c.cent9Max + 1; ++iz) {
        sum += h->GetBinContent(ix, iy, iz);
        const Double_t e = h->GetBinError(ix, iy, iz);
        err2 += e * e;
      }
      out->SetBinContent(ix, iy, sum);
      out->SetBinError(ix, iy, TMath::Sqrt(err2));
    }
  }
  return out;
}

void MeanMtByKstar(TH3F* h, const Data006CentralityConfig& c, Double_t kmax,
                   const std::string& channel, const std::string& sample,
                   std::ofstream& csv) {
  if (!h) return;
  const Int_t nk = h->GetYaxis()->FindFixBin(kmax - 1e-12);
  for (Int_t ik = 1; ik <= nk; ++ik) {
    Double_t n = 0.0, s = 0.0, s2 = 0.0;
    for (Int_t im = 0; im <= h->GetXaxis()->GetNbins() + 1; ++im) {
      Double_t w = 0.0;
      for (Int_t iz = c.cent9Min + 1; iz <= c.cent9Max + 1; ++iz)
        w += h->GetBinContent(im, ik, iz);
      const Double_t x = h->GetXaxis()->GetBinCenter(im);
      n += w; s += w * x; s2 += w * x * x;
    }
    const Double_t mean = n > 0.0 ? s / n : QuietNaN();
    const Double_t var = n > 0.0 ? TMath::Max(0.0, s2 / n - mean * mean) : QuietNaN();
    const Double_t rms = n > 0.0 ? TMath::Sqrt(var) : QuietNaN();
    const Double_t sem = n > 0.0 ? rms / TMath::Sqrt(n) : QuietNaN();
    csv << channel << ',' << c.id << ',' << sample << ','
        << 1000.0 * h->GetYaxis()->GetBinLowEdge(ik) << ','
        << 1000.0 * h->GetYaxis()->GetBinUpEdge(ik) << ','
        << n << ',' << mean << ',' << rms << ',' << sem << ','
        << (n > 0.0 ? "OK" : "EMPTY") << '\n';
  }
}

}  // namespace

void finalizeData006(const Char_t* mergedInput, const Char_t* data006Config,
                     const Char_t* outputDirectory, const Char_t* productionTag,
                     const Char_t* gitCommit) {
  Data006Config cfg;
  if (!cfg.Load(data006Config)) {
    std::cerr << "ERROR: invalid DATA-006 config: " << data006Config << std::endl;
    return;
  }
  TFile* in = TFile::Open(mergedInput, "READ");
  if (!in || in->IsZombie()) {
    std::cerr << "ERROR: cannot open merged DATA-006 input " << mergedInput << std::endl;
    return;
  }
  if (gSystem->mkdir(outputDirectory, kTRUE) != 0 && gSystem->AccessPathName(outputDirectory)) {
    std::cerr << "ERROR: cannot create output directory " << outputDirectory << std::endl;
    return;
  }

  TH1D* hEvents = (TH1D*)in->Get("hAcceptedEventsVsCent9");
  TH1D* hMult = (TH1D*)in->Get("hFxtMultSumVsCent9");
  TH1D* hMult2 = (TH1D*)in->Get("hFxtMultSum2VsCent9");
  if (!hEvents || !hMult || !hMult2) {
    std::cerr << "ERROR: merged input lacks exact DATA-006 event/multiplicity accumulators"
              << std::endl;
    return;
  }

  const TString rootPath = TString::Format("%s/data006_correlations.root", outputDirectory);
  const TString corrPath = TString::Format("%s/correlations.csv", outputDirectory);
  const TString centPath = TString::Format("%s/centrality_summary.csv", outputDirectory);
  const TString mtPath = TString::Format("%s/pair_mt_summary.csv", outputDirectory);
  const TString mtKPath = TString::Format("%s/pair_mt_by_kstar.csv", outputDirectory);
  const TString manifestPath = TString::Format("%s/object_manifest.csv", outputDirectory);
  const TString invalidPath = TString::Format("%s/invalid_bins.csv", outputDirectory);
  const TString qaPath = TString::Format("%s/data006_correlation_qa.pdf", outputDirectory);

  TFile* out = TFile::Open(rootPath, "RECREATE");
  if (!out || out->IsZombie()) {
    std::cerr << "ERROR: cannot create " << rootPath << std::endl;
    return;
  }
  std::ofstream corr(corrPath.Data()), cent(centPath.Data()), mt(mtPath.Data());
  std::ofstream mtK(mtKPath.Data()), manifest(manifestPath.Data()), invalid(invalidPath.Data());
  if (!corr || !cent || !mt || !mtK || !manifest || !invalid) {
    std::cerr << "ERROR: cannot create one or more DATA-006 CSV outputs" << std::endl;
    return;
  }
  corr << std::setprecision(17);
  cent << std::setprecision(17);
  mt << std::setprecision(17);
  mtK << std::setprecision(17);
  invalid << std::setprecision(17);

  corr << "channel,centrality,bin_kind,kstar_low_mev_c,kstar_high_mev_c,se,se_error,me,me_error,normalization_factor,normalization_error,normalization_min_mev_c,normalization_max_mev_c,cf,cf_error,status\n";
  cent << "centrality,low_percent,high_percent,cent9_min,cent9_max,native_class,refmultcorr_low_exclusive,refmultcorr_high_inclusive,event_count,fxtmult_raw_mean,fxtmult_raw_rms,fxtmult_raw_sem,event_weighting,nch_corr_mean,nch_corr_rms,nch_syst_error,production_tag,git_commit\n";
  mt << "channel,centrality,sample,pair_count,mean_pair_mt_half_gev_c2,rms_pair_mt_half_gev_c2,sem_pair_mt_half_gev_c2,kstar_max_mev_c,status\n";
  mtK << "channel,centrality,sample,kstar_low_mev_c,kstar_high_mev_c,pair_count,mean_pair_mt_half_gev_c2,rms_pair_mt_half_gev_c2,sem_pair_mt_half_gev_c2,status\n";
  manifest << "channel,centrality,observable,root_file,root_object,csv_file,unit,normalization,binning\n";
  invalid << "channel,centrality,bin_kind,kstar_low_mev_c,kstar_high_mev_c,reason\n";

  Bool_t allRequired = kTRUE;
  Long64_t invalidBins = 0;
  for (size_t iz = 0; iz < cfg.centralities.size(); ++iz) {
    const Data006CentralityConfig& z = cfg.centralities[iz];
    Double_t n = 0.0, s = 0.0, s2 = 0.0;
    for (Int_t ic = z.cent9Min + 1; ic <= z.cent9Max + 1; ++ic) {
      n += hEvents->GetBinContent(ic);
      s += hMult->GetBinContent(ic);
      s2 += hMult2->GetBinContent(ic);
    }
    const Double_t mean = n > 0.0 ? s / n : QuietNaN();
    const Double_t var = n > 0.0 ? TMath::Max(0.0, s2 / n - mean * mean) : QuietNaN();
    const Double_t rms = n > 0.0 ? TMath::Sqrt(var) : QuietNaN();
    const Double_t sem = n > 0.0 ? rms / TMath::Sqrt(n) : QuietNaN();
    cent << z.id << ',' << z.lowPercent << ',' << z.highPercent << ',' << z.cent9Min << ','
         << z.cent9Max << ',' << (z.nativeClass ? "true" : "false") << ',' << z.refmultLow
         << ',' << z.refmultHigh << ',' << n << ',' << mean << ',' << rms << ',' << sem
         << ",unweighted_accepted_events,NOT_AVAILABLE,NOT_AVAILABLE,NOT_AVAILABLE,\""
         << productionTag << "\"," << gitCommit << '\n';
  }

  Bool_t firstQa = kTRUE;
  for (size_t ich = 0; ich < cfg.channels.size(); ++ich) {
    const Data006ChannelConfig& ch = cfg.channels[ich];
    if (!ch.enabled) continue;
    TH2D* seCent = (TH2D*)in->Get(TString::Format("hKstarSEVsCent_%s", ch.name.c_str()));
    TH2D* meCent = (TH2D*)in->Get(TString::Format("hKstarMEVsCent_%s", ch.name.c_str()));
    TH3F* mtSE = (TH3F*)in->Get(TString::Format("hPairMtHalf_vs_KstarSEVsCent_%s", ch.name.c_str()));
    TH3F* mtME = (TH3F*)in->Get(TString::Format("hPairMtHalf_vs_KstarMEVsCent_%s", ch.name.c_str()));
    TH2D* ySE = (TH2D*)in->Get(TString::Format("hPairRapiditySEVsCent_%s", ch.name.c_str()));
    TH2D* yME = (TH2D*)in->Get(TString::Format("hPairRapidityMEVsCent_%s", ch.name.c_str()));
    TH2D* yNormSE = (TH2D*)in->Get(
        TString::Format("hPairRapidityNormSEVsCent_%s", ch.name.c_str()));
    TH2D* yNormME = (TH2D*)in->Get(
        TString::Format("hPairRapidityNormMEVsCent_%s", ch.name.c_str()));
    TH3F* closeSE = (TH3F*)in->Get(
        TString::Format("hDeltaEtaDeltaPhiStarSEVsCent_%s", ch.name.c_str()));
    TH3F* closeME = (TH3F*)in->Get(
        TString::Format("hDeltaEtaDeltaPhiStarMEVsCent_%s", ch.name.c_str()));
    TH1D* cutFlow = (TH1D*)in->Get(TString::Format("hPairCutFlow_%s", ch.name.c_str()));
    const Char_t* kinName[8] = {"PartAPt","PartAP","PartARapidity","PartAEta",
                                "PartBPt","PartBP","PartBRapidity","PartBEta"};
    TH2D* kinSourceSE[8] = {0};
    TH2D* kinSourceME[8] = {0};
    for (Int_t ik = 0; ik < 8; ++ik) {
      kinSourceSE[ik] = (TH2D*)in->Get(
          TString::Format("h%sSEVsCent_%s", kinName[ik], ch.name.c_str()));
      kinSourceME[ik] = (TH2D*)in->Get(
          TString::Format("h%sMEVsCent_%s", kinName[ik], ch.name.c_str()));
    }
    Bool_t channelObjectsOk = seCent && meCent && mtSE && mtME && ySE && yME &&
                              yNormSE && yNormME && closeSE && closeME && cutFlow;
    for (Int_t ik = 0; ik < 8; ++ik)
      channelObjectsOk = channelObjectsOk && kinSourceSE[ik] && kinSourceME[ik];
    if (!channelObjectsOk) {
      std::cerr << "ERROR: missing required merged objects for " << ch.name << std::endl;
      allRequired = kFALSE;
      continue;
    }

    TCanvas qa(TString::Format("qa_%s", ch.name.c_str()), ch.name.c_str(), 1800, 900);
    TCanvas qaAcc(TString::Format("qaAcc_%s", ch.name.c_str()), ch.name.c_str(), 1800, 900);
    qaAcc.Divide(5, 2);
    qa.Divide(5, 2);
    Int_t pad = 1;
    Int_t accPad = 1;
    for (size_t iz = 0; iz < cfg.centralities.size(); ++iz) {
      const Data006CentralityConfig& z = cfg.centralities[iz];
      TString base = TString::Format("%s_%s", ch.name.c_str(), z.id.c_str());
      TH1D* se = ProjectCent(seCent, z, TString::Format("hSE_%s", base.Data()));
      TH1D* me = ProjectCent(meCent, z, TString::Format("hME_%s", base.Data()));
      if (!se || !me) { allRequired = kFALSE; continue; }
      const Int_t n1 = se->GetXaxis()->FindFixBin(ch.normQMin + 1e-12);
      const Int_t n2 = se->GetXaxis()->FindFixBin(ch.normQMax - 1e-12);
      Double_t seNormErr = 0.0, meNormErr = 0.0;
      const Double_t seNorm = se->IntegralAndError(n1, n2, seNormErr);
      const Double_t meNorm = me->IntegralAndError(n1, n2, meNormErr);
      const Double_t alpha = (meNorm > 0.0) ? seNorm / meNorm : QuietNaN();
      const Double_t alphaErr = (seNorm > 0.0 && meNorm > 0.0)
          ? alpha * TMath::Sqrt((seNormErr / seNorm) * (seNormErr / seNorm) +
                               (meNormErr / meNorm) * (meNormErr / meNorm))
          : QuietNaN();
      if (!TMath::Finite(alpha) || alpha <= 0.0) allRequired = kFALSE;

      TH1D* cf = (TH1D*)se->Clone(TString::Format("hCF_%s", base.Data()));
      cf->Reset();
      cf->SetTitle(TString::Format("%s %d-%d%%;k* (GeV/c);C(k*)", ch.name.c_str(),
                                   z.lowPercent, z.highPercent));
      for (Int_t ib = 0; ib <= se->GetNbinsX() + 1; ++ib) {
        const Double_t sv = se->GetBinContent(ib), sev = se->GetBinError(ib);
        const Double_t mv = me->GetBinContent(ib), mev = me->GetBinError(ib);
        TString kind("native");
        if (ib == 0) kind = "underflow";
        else if (ib == se->GetNbinsX() + 1) kind = "overflow";
        const Double_t lo = 1000.0 * se->GetXaxis()->GetBinLowEdge(ib);
        const Double_t hi = 1000.0 * se->GetXaxis()->GetBinUpEdge(ib);
        Double_t cv = QuietNaN(), ce = QuietNaN();
        TString status("OK");
        const Bool_t deliver = ib == 0 || ib == se->GetNbinsX() + 1 ||
                               se->GetXaxis()->GetBinLowEdge(ib) < cfg.deliveryKstarMax;
        if (mv <= 0.0 || !TMath::Finite(alpha) || alpha <= 0.0) {
          status = mv <= 0.0 ? "ZERO_DENOMINATOR" : "INVALID_NORMALIZATION";
          if (deliver) {
            ++invalidBins;
            invalid << ch.name << ',' << z.id << ',' << kind << ',' << lo << ',' << hi << ','
                    << status << '\n';
          }
        } else {
          cv = sv / (alpha * mv);
          const Double_t dS = 1.0 / (alpha * mv);
          const Double_t dM = -sv / (alpha * mv * mv);
          const Double_t dA = -sv / (alpha * alpha * mv);
          ce = TMath::Sqrt(dS * dS * sev * sev + dM * dM * mev * mev +
                           dA * dA * alphaErr * alphaErr);
        }
        cf->SetBinContent(ib, cv);
        cf->SetBinError(ib, ce);
        if (deliver)
          corr << ch.name << ',' << z.id << ',' << kind << ',' << lo << ',' << hi << ','
               << sv << ',' << sev << ',' << mv << ',' << mev << ',' << alpha << ','
               << alphaErr << ',' << 1000.0 * ch.normQMin << ',' << 1000.0 * ch.normQMax
               << ',' << cv << ',' << ce << ',' << status << '\n';
      }

      TH1D* mtSe1 = ProjectMt(mtSE, z, cfg.deliveryKstarMax,
                              TString::Format("hPairMtHalfSE_%s", base.Data()));
      TH1D* mtMe1 = ProjectMt(mtME, z, cfg.deliveryKstarMax,
                              TString::Format("hPairMtHalfME_%s", base.Data()));
      TH2D* mtSe2 = ProjectKstarMt(mtSE, z, cfg.deliveryKstarMax,
                                   TString::Format("hKstarVsPairMtHalfSE_%s", base.Data()));
      TH2D* mtMe2 = ProjectKstarMt(mtME, z, cfg.deliveryKstarMax,
                                   TString::Format("hKstarVsPairMtHalfME_%s", base.Data()));
      TH1D* ySe1 = ProjectPairY(ySE, z, TString::Format("hPairRapiditySE_%s", base.Data()));
      TH1D* yMe1 = ProjectPairY(yME, z, TString::Format("hPairRapidityME_%s", base.Data()));
      TH1D* yNormSe1 = ProjectPairY(
          yNormSE, z, TString::Format("hPairRapidityNormSE_%s", base.Data()));
      TH1D* yNormMe1 = ProjectPairY(
          yNormME, z, TString::Format("hPairRapidityNormME_%s", base.Data()));
      TH2D* closeSe2 = ProjectXYCent(
          closeSE, z, TString::Format("hClosePairQA_SE_%s", base.Data()),
          "SE two-track QA;#Delta#eta;min #Delta#phi* (rad)");
      TH2D* closeMe2 = ProjectXYCent(
          closeME, z, TString::Format("hClosePairQA_ME_%s", base.Data()),
          "ME two-track QA;#Delta#eta;min #Delta#phi* (rad)");
      TH1D* kinSe1[8] = {0};
      TH1D* kinMe1[8] = {0};
      for (Int_t ik = 0; ik < 8; ++ik) {
        kinSe1[ik] = ProjectPairY(
            kinSourceSE[ik], z, TString::Format("h%sSE_%s", kinName[ik], base.Data()));
        kinMe1[ik] = ProjectPairY(
            kinSourceME[ik], z, TString::Format("h%sME_%s", kinName[ik], base.Data()));
      }

      TDirectory* dch = out->GetDirectory(ch.name.c_str());
      if (!dch) dch = out->mkdir(ch.name.c_str());
      TDirectory* dz = dch->GetDirectory(z.id.c_str());
      if (!dz) dz = dch->mkdir(z.id.c_str());
      dz->cd();
      se->Write("hSE"); me->Write("hME"); cf->Write("hCF");
      mtSe1->Write("hPairMtHalfSE"); mtMe1->Write("hPairMtHalfME");
      mtSe2->Write("hKstarVsPairMtHalfSE"); mtMe2->Write("hKstarVsPairMtHalfME");
      ySe1->Write("hPairRapiditySE"); yMe1->Write("hPairRapidityME");
      yNormSe1->Write("hPairRapidityNormSE");
      yNormMe1->Write("hPairRapidityNormME");
      closeSe2->Write("hClosePairQA_SE");
      closeMe2->Write("hClosePairQA_ME");
      cutFlow->Write("hPairCutFlow");
      for (Int_t ik = 0; ik < 8; ++ik)
        kinSe1[ik]->Write(TString::Format("h%sSE", kinName[ik]));
      for (Int_t ik = 0; ik < 8; ++ik)
        kinMe1[ik]->Write(TString::Format("h%sME", kinName[ik]));
      TParameter<Double_t>("normalizationFactor", alpha).Write();
      TParameter<Double_t>("normalizationError", alphaErr).Write();
      TParameter<Double_t>("normalizationMinGeV", ch.normQMin).Write();
      TParameter<Double_t>("normalizationMaxGeV", ch.normQMax).Write();
      TNamed("mergedClassConstruction", "SE and ME projected/summed first; normalization and CF recomputed afterward").Write();
      out->cd();

      mt << ch.name << ',' << z.id << ",SE," << mtSe1->GetEntries() << ','
         << mtSe1->GetMean() << ',' << mtSe1->GetRMS() << ',' << mtSe1->GetMeanError() << ','
         << 1000.0 * cfg.deliveryKstarMax << ',' << (mtSe1->GetEntries() > 0 ? "OK" : "EMPTY") << '\n';
      mt << ch.name << ',' << z.id << ",ME," << mtMe1->GetEntries() << ','
         << mtMe1->GetMean() << ',' << mtMe1->GetRMS() << ',' << mtMe1->GetMeanError() << ','
         << 1000.0 * cfg.deliveryKstarMax << ',' << (mtMe1->GetEntries() > 0 ? "OK" : "EMPTY") << '\n';
      MeanMtByKstar(mtSE, z, cfg.deliveryKstarMax, ch.name, "SE", mtK);
      MeanMtByKstar(mtME, z, cfg.deliveryKstarMax, ch.name, "ME", mtK);

      const Double_t kWidth = (cfg.kstarMax - cfg.kstarMin) / cfg.kstarBins;
      const Int_t deliveredKBins = (Int_t)((cfg.deliveryKstarMax - cfg.kstarMin) / kWidth + 0.5);
      const Char_t* samples[2] = {"SE", "ME"};
      for (Int_t is = 0; is < 2; ++is)
        manifest << ch.name << ',' << z.id << ',' << samples[is]
                 << ",data006_correlations.root," << ch.name << '/' << z.id << "/h"
                 << samples[is] << ",correlations.csv,pairs,kstar_raw,kstar_GeV_c:"
                 << cfg.kstarMin << ':' << cfg.kstarMax << ':' << cfg.kstarBins << '\n';
      manifest << ch.name << ',' << z.id << ",CF,data006_correlations.root," << ch.name
               << '/' << z.id
               << "/hCF,correlations.csv,dimensionless,SE/(alpha*ME),kstar_GeV_c:"
               << cfg.kstarMin << ':' << cfg.kstarMax << ':' << cfg.kstarBins << '\n';
      for (Int_t is = 0; is < 2; ++is)
        manifest << ch.name << ',' << z.id << ",pair_mt_half_"
                 << (is == 0 ? "se" : "me") << ",data006_correlations.root," << ch.name
                 << '/' << z.id << "/hPairMtHalf" << samples[is]
                 << ",pair_mt_summary.csv,GeV/c^2,none,pair_mt_half_GeV_c2:"
                 << cfg.pairMtMin << ':' << cfg.pairMtMax << ':' << cfg.pairMtBins << '\n';
      for (Int_t is = 0; is < 2; ++is)
        manifest << ch.name << ',' << z.id << ",kstar_vs_pair_mt_half_"
                 << (is == 0 ? "se" : "me") << ",data006_correlations.root," << ch.name
                 << '/' << z.id << "/hKstarVsPairMtHalf" << samples[is]
                 << ",pair_mt_by_kstar.csv,GeV/c_and_GeV/c^2,none,kstar_GeV_c:0:"
                 << cfg.deliveryKstarMax << ':' << deliveredKBins << "|pair_mt_half_GeV_c2:"
                 << cfg.pairMtMin << ':' << cfg.pairMtMax << ':' << cfg.pairMtBins << '\n';
      for (Int_t is = 0; is < 2; ++is) {
        manifest << ch.name << ',' << z.id << ",pair_rapidity_"
                 << (is == 0 ? "se" : "me") << ",data006_correlations.root," << ch.name
                 << '/' << z.id << "/hPairRapidity" << samples[is]
                 << ",NOT_AVAILABLE,dimensionless,none,pair_y_cm:"
                 << cfg.pairRapidityMin << ':' << cfg.pairRapidityMax << ':'
                 << cfg.pairRapidityBins << '\n';
        manifest << ch.name << ',' << z.id << ",pair_rapidity_norm_"
                 << (is == 0 ? "se" : "me") << ",data006_correlations.root," << ch.name
                 << '/' << z.id << "/hPairRapidityNorm" << samples[is]
                 << ",NOT_AVAILABLE,dimensionless,norm_kstar_only,pair_y_cm:"
                 << cfg.pairRapidityMin << ':' << cfg.pairRapidityMax << ':'
                 << cfg.pairRapidityBins << '\n';
        manifest << ch.name << ',' << z.id << ",close_pair_qa_"
                 << (is == 0 ? "se" : "me") << ",data006_correlations.root," << ch.name
                 << '/' << z.id << "/hClosePairQA_" << samples[is]
                 << ",NOT_AVAILABLE,dimensionless_and_rad,none,delta_eta:-0.16:0.16:160|delta_phi_star_rad:-0.16:0.16:160\n";
      }
      manifest << ch.name << ',' << z.id
               << ",pair_cut_flow,data006_correlations.root," << ch.name << '/' << z.id
               << "/hPairCutFlow,NOT_AVAILABLE,pairs,none,categorical:10\n";
      for (Int_t ik = 0; ik < 8; ++ik) {
        TString kinBins;
        if (ik == 0 || ik == 1 || ik == 4 || ik == 5) kinBins = "axis_GeV_c:0:5:250";
        else if (ik == 2 || ik == 6)
          kinBins = TString::Format("axis_y_cm:%.12g:%.12g:%d", cfg.pairRapidityMin,
                                    cfg.pairRapidityMax, cfg.pairRapidityBins);
        else kinBins = "axis_eta:-2.5:2.5:200";
        for (Int_t is = 0; is < 2; ++is)
          manifest << ch.name << ',' << z.id << ',' << kinName[ik]
                   << (is == 0 ? "_se" : "_me") << ",data006_correlations.root," << ch.name
                   << '/' << z.id << "/h" << kinName[ik] << samples[is]
                   << ",NOT_AVAILABLE,configured_axis,none," << kinBins << '\n';
      }
      manifest << ch.name << ',' << z.id
               << ",normalization_factor,data006_correlations.root," << ch.name << '/' << z.id
               << "/normalizationFactor,correlations.csv,dimensionless,Integral(SE)/Integral(ME),scalar\n";

      if (accPad <= 10) {
        qaAcc.cd(accPad++);
        TH1D* yNormMeDraw = (TH1D*)yNormMe1->Clone(
            TString::Format("hPairRapidityNormMEDraw_%s", base.Data()));
        const Double_t seY = yNormSe1->Integral();
        const Double_t meY = yNormMeDraw->Integral();
        if (seY > 0.0 && meY > 0.0) yNormMeDraw->Scale(seY / meY);
        yNormSe1->SetLineColor(1);
        yNormMeDraw->SetLineColor(2);
        yNormSe1->SetTitle(TString::Format(
            "%s %d-%d%% norm-k* pair y: SE black, scaled ME red;y_{CM};pairs",
            ch.name.c_str(), z.lowPercent, z.highPercent));
        yNormSe1->Draw("HIST");
        yNormMeDraw->Draw("HIST SAME");
      }
      if (pad <= 10) {
        qa.cd(pad++);
        cf->GetXaxis()->SetRangeUser(0.0, cfg.deliveryKstarMax);
        cf->SetMinimum(0.0); cf->SetMaximum(2.5); cf->Draw("E1");
      }
    }
    TString cfPrint = qaPath;
    if (firstQa) cfPrint += "(";
    qa.Print(cfPrint);
    Bool_t lastQa = kTRUE;
    for (size_t jch = ich + 1; jch < cfg.channels.size(); ++jch)
      if (cfg.channels[jch].enabled) lastQa = kFALSE;
    TString accPrint = qaPath;
    if (lastQa) accPrint += ")";
    qaAcc.Print(accPrint);
    firstQa = kFALSE;
  }

  out->cd();
  TNamed("inputMergedRoot", mergedInput).Write();
  TNamed("data006Config", data006Config).Write();
  TNamed("productionTag", productionTag).Write();
  TNamed("gitCommit", gitCommit).Write();
  TNamed("kstarUnit", "GeV/c in ROOT; MeV/c in CSV").Write();
  TNamed("normalizationMethod", "alpha=Integral(SE)/Integral(ME); C=SE/(alpha*ME)").Write();
  TNamed("normalizationErrorMethod", "independent integral errors propagated to alpha and each CF bin; norm-bin covariance not subtracted").Write();
  TNamed("nchCorr", "NOT_AVAILABLE").Write();
  TNamed("purityPrimaryFeeddownSystematics", "NOT_AVAILABLE").Write();
  TParameter<Long64_t>("acceptedEventCount", (Long64_t)hEvents->Integral(1, 9)).Write();
  TParameter<Long64_t>("invalidBinCount", invalidBins).Write();
  TNamed("closureStatus", allRequired ? "FIT_READY" : "NOT_FIT_READY").Write();
  out->Close();
  in->Close();

  std::ofstream status(TString::Format("%s/CLOSURE_STATUS.txt", outputDirectory).Data());
  status << (allRequired ? "FIT_READY" : "NOT_FIT_READY") << '\n'
         << "invalid_or_zero_denominator_bins=" << invalidBins << '\n'
         << "missing systematic inputs are enumerated in ROOT metadata and README.\n";
  std::cout << "[finalizeData006] " << (allRequired ? "FIT_READY" : "NOT_FIT_READY")
            << " output=" << outputDirectory << " invalidBins=" << invalidBins << std::endl;
}
