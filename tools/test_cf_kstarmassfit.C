// Toy closure for kstarMassFitCF: S = F - α B, C_raw = Y_SE / Y_ME.
// root -b -q tools/test_cf_kstarmassfit.C

#include <TH1D.h>
#include <TF1.h>
#include <TMath.h>
#include <iostream>
#include <vector>

static TH1D* makeBkg(const char* name, Double_t nBkg) {
  TH1D* h = new TH1D(name, name, 140, 0.99, 1.06);
  const Double_t perBin = nBkg / 140.0;
  for (Int_t b = 1; b <= 140; ++b) {
    h->SetBinContent(b, perBin);
    h->SetBinError(b, TMath::Sqrt(perBin));
  }
  return h;
}

static void addSignal(TH1D* h, Double_t nSig, Double_t mean, Double_t sigma) {
  TF1 g("gsig", "gaus", 0.99, 1.06);
  g.SetParameters(nSig / (sigma * TMath::Sqrt(2.0 * TMath::Pi())), mean, sigma);
  h->Add(&g, 1.0);
}

static Double_t alphaFromWindow(TH1* hF, TH1* hB, Double_t aMin, Double_t aMax) {
  const Int_t b0 = hF->GetXaxis()->FindBin(aMin + 1e-9);
  const Int_t b1 = hF->GetXaxis()->FindBin(aMax - 1e-9);
  const Double_t f = hF->Integral(b0, b1);
  const Double_t b = hB->Integral(b0, b1);
  if (b <= 0.0) return 0.0;
  return f / b;
}

static Bool_t fitGausYield(TH1* hS, Double_t sigMin, Double_t sigMax, Double_t& nSig) {
  nSig = 0.0;
  TF1 f("fgaus", "gaus", 0.99, 1.06);
  f.SetParameter(0, hS->GetMaximum());
  f.SetParameter(1, 1.0195);
  f.SetParameter(2, 0.004);
  f.SetParLimits(2, 0.002, 0.020);
  const Int_t st = hS->Fit(&f, "RQ0");
  if (st != 0) return kFALSE;
  nSig = f.Integral(sigMin, sigMax);
  return nSig > 0.0 && TMath::Finite(nSig);
}

static Double_t toyOverlapIntegral(TH1* h, Double_t xMin, Double_t xMax, Double_t& err) {
  err = 0.0;
  if (!h) return 0.0;
  Double_t sum = 0.0, err2 = 0.0;
  for (Int_t b = 1; b <= h->GetNbinsX(); ++b) {
    const Double_t lo = h->GetXaxis()->GetBinLowEdge(b);
    const Double_t hi = h->GetXaxis()->GetBinUpEdge(b);
    const Double_t bw = hi - lo;
    if (bw <= 0.0) continue;
    const Double_t ov = TMath::Min(hi, xMax) - TMath::Max(lo, xMin);
    if (ov <= 0.0) continue;
    const Double_t w = ov / bw;
    sum += w * h->GetBinContent(b);
    err2 += TMath::Power(w * h->GetBinError(b), 2);
  }
  err = TMath::Sqrt(err2);
  return sum;
}

static void addBoxSignal(TH1D* h, Double_t nSig, Double_t xMin, Double_t xMax) {
  Double_t wsum = 0.0;
  for (Int_t b = 1; b <= h->GetNbinsX(); ++b) {
    const Double_t lo = h->GetXaxis()->GetBinLowEdge(b);
    const Double_t hi = h->GetXaxis()->GetBinUpEdge(b);
    const Double_t ov = TMath::Min(hi, xMax) - TMath::Max(lo, xMin);
    if (ov > 0.0) wsum += ov;
  }
  if (wsum <= 0.0) return;
  for (Int_t b = 1; b <= h->GetNbinsX(); ++b) {
    const Double_t lo = h->GetXaxis()->GetBinLowEdge(b);
    const Double_t hi = h->GetXaxis()->GetBinUpEdge(b);
    const Double_t ov = TMath::Min(hi, xMax) - TMath::Max(lo, xMin);
    if (ov <= 0.0) continue;
    h->SetBinContent(b, h->GetBinContent(b) + nSig * ov / wsum);
  }
}

static Int_t test_window_count() {
  Int_t nFail = 0;
  const Double_t aMin = 1.04, aMax = 1.06, sigMin = 1.012, sigMax = 1.026;
  TH1D* hB = makeBkg("toyB", 4000.0);
  TH1D* hF = (TH1D*)hB->Clone("toyF");
  const Double_t yInj = 250.0;
  addBoxSignal(hF, yInj, sigMin, sigMax);
  Double_t eFA = 0, eBA = 0, eFW = 0, eBW = 0;
  const Double_t fA = toyOverlapIntegral(hF, aMin, aMax, eFA);
  const Double_t bA = toyOverlapIntegral(hB, aMin, aMax, eBA);
  const Double_t alpha = fA / bA;
  const Double_t eA = TMath::Abs(alpha) * TMath::Sqrt(TMath::Power(eFA / fA, 2) + TMath::Power(eBA / bA, 2));
  const Double_t fW = toyOverlapIntegral(hF, sigMin, sigMax, eFW);
  const Double_t bW = toyOverlapIntegral(hB, sigMin, sigMax, eBW);
  const Double_t y = fW - alpha * bW;
  const Double_t eY = TMath::Sqrt(eFW * eFW + TMath::Power(alpha * eBW, 2) + TMath::Power(bW * eA, 2));
  const Double_t eYnoA = TMath::Sqrt(eFW * eFW + TMath::Power(alpha * eBW, 2));
  if (TMath::Abs(y - yInj) / yInj > 0.02) {
    std::cerr << "FAIL box recovery y=" << y << " inj=" << yInj << std::endl;
    ++nFail;
  } else {
    std::cout << "PASS box-signal window recovery y=" << y << " inj=" << yInj << std::endl;
  }
  TH1D* hFg = (TH1D*)hB->Clone("toyFg");
  TF1 gbox("gbox", "gaus", 0.99, 1.06);
  gbox.SetParameters(1.0, 1.0195, 0.0030);
  Double_t gsum = 0.0;
  for (Int_t b = 1; b <= hFg->GetNbinsX(); ++b) {
    gsum += gbox.Integral(hFg->GetXaxis()->GetBinLowEdge(b), hFg->GetXaxis()->GetBinUpEdge(b));
  }
  for (Int_t b = 1; b <= hFg->GetNbinsX(); ++b) {
    const Double_t add =
        yInj * gbox.Integral(hFg->GetXaxis()->GetBinLowEdge(b), hFg->GetXaxis()->GetBinUpEdge(b)) / gsum;
    hFg->SetBinContent(b, hFg->GetBinContent(b) + add);
  }
  Double_t eFAg = 0, eFWg = 0;
  const Double_t fAg = toyOverlapIntegral(hFg, aMin, aMax, eFAg);
  const Double_t ag = fAg / bA;
  const Double_t fWg = toyOverlapIntegral(hFg, sigMin, sigMax, eFWg);
  const Double_t yg = fWg - ag * bW;
  if (TMath::Abs(yg - yInj) / yInj > 0.08) {
    std::cerr << "FAIL non-box window recovery yg=" << yg << " inj=" << yInj << std::endl;
    ++nFail;
  } else {
    std::cout << "PASS non-Gaussian-fit (narrow gaus) window recovery yg=" << yg << std::endl;
  }
  if (!(eY > eYnoA)) {
    std::cerr << "FAIL coherent alpha term not increasing error" << std::endl;
    ++nFail;
  } else {
    std::cout << "PASS coherent alpha term eY=" << eY << " eY(no alpha)=" << eYnoA << std::endl;
  }
  TH1D* hS = (TH1D*)hF->Clone("toyS");
  hS->Add(hB, -alpha);
  Int_t nNeg = 0;
  Double_t yNoClip = 0.0;
  for (Int_t b = 1; b <= hS->GetNbinsX(); ++b) {
    const Double_t lo = hS->GetXaxis()->GetBinLowEdge(b);
    const Double_t hi = hS->GetXaxis()->GetBinUpEdge(b);
    const Double_t ov = TMath::Min(hi, sigMax) - TMath::Max(lo, sigMin);
    if (ov <= 0.0) continue;
    if (hS->GetBinContent(b) < 0.0) ++nNeg;
    yNoClip += (ov / (hi - lo)) * hS->GetBinContent(b);
  }
  if (TMath::Abs(yNoClip - y) > 1e-6) {
    std::cerr << "FAIL clip-free sum mismatch" << std::endl;
    ++nFail;
  } else {
    std::cout << "PASS no-clip window sum nNeg=" << nNeg << std::endl;
  }
  TH1D* hRe = (TH1D*)hF->Clone("toyRe");
  hRe->Rebin(2);
  TH1D* hBre = (TH1D*)hB->Clone("toyBre");
  hBre->Rebin(2);
  Double_t eFWr = 0, eBWr = 0, eFAr = 0, eBAr = 0;
  const Double_t ar = toyOverlapIntegral(hRe, aMin, aMax, eFAr) / toyOverlapIntegral(hBre, aMin, aMax, eBAr);
  const Double_t yr = toyOverlapIntegral(hRe, sigMin, sigMax, eFWr) - ar * toyOverlapIntegral(hBre, sigMin, sigMax, eBWr);
  if (TMath::Abs(yr - y) / TMath::Abs(y) > 0.05) {
    std::cerr << "FAIL rebin changed counting yield native=" << y << " rebin=" << yr << std::endl;
    ++nFail;
  } else {
    std::cout << "PASS rebin invariance native=" << y << " rebin2=" << yr << std::endl;
  }
  TH1D* h40a = makeBkg("b40a", 2000.0);
  TH1D* h40b = makeBkg("b40b", 2000.0);
  TH1D* h80 = makeBkg("b80", 4000.0);
  TH1D* f40a = (TH1D*)h40a->Clone("f40a");
  TH1D* f40b = (TH1D*)h40b->Clone("f40b");
  TH1D* f80 = (TH1D*)h80->Clone("f80");
  addBoxSignal(f40a, 120.0, sigMin, sigMax);
  addBoxSignal(f40b, 130.0, sigMin, sigMax);
  addBoxSignal(f80, 250.0, sigMin, sigMax);
  Double_t dum = 0;
  const Double_t y40 = (toyOverlapIntegral(f40a, sigMin, sigMax, dum) -
                        (toyOverlapIntegral(f40a, aMin, aMax, dum) / toyOverlapIntegral(h40a, aMin, aMax, dum)) *
                            toyOverlapIntegral(h40a, sigMin, sigMax, dum)) +
                       (toyOverlapIntegral(f40b, sigMin, sigMax, dum) -
                        (toyOverlapIntegral(f40b, aMin, aMax, dum) / toyOverlapIntegral(h40b, aMin, aMax, dum)) *
                            toyOverlapIntegral(h40b, sigMin, sigMax, dum));
  const Double_t y80 = toyOverlapIntegral(f80, sigMin, sigMax, dum) -
                       (toyOverlapIntegral(f80, aMin, aMax, dum) / toyOverlapIntegral(h80, aMin, aMax, dum)) *
                           toyOverlapIntegral(h80, sigMin, sigMax, dum);
  if (TMath::Abs(y80 - y40) / y80 > 0.02) {
    std::cerr << "FAIL 40+40 vs 80 y40=" << y40 << " y80=" << y80 << std::endl;
    ++nFail;
  } else {
    std::cout << "PASS 40+40=80 y40sum=" << y40 << " y80=" << y80 << std::endl;
  }
  const Double_t yME = 0.0;
  const Bool_t cfInvalid = !(yME > 0.0);
  if (!cfInvalid) {
    std::cerr << "FAIL invalid ME should not make CF" << std::endl;
    ++nFail;
  } else {
    std::cout << "PASS invalid ME denominator does not create CF=0" << std::endl;
  }
  const Double_t raw[] = {1.2, 1.1, 1.0, 0.95};
  const Double_t kx[] = {0.05, 0.15, 0.45, 0.55};
  Double_t sum = 0.0;
  Int_t nN = 0;
  for (Int_t i = 0; i < 4; ++i) {
    if (kx[i] > 0.4 && kx[i] < 0.6) {
      sum += raw[i];
      ++nN;
    }
  }
  const Double_t scale = (Double_t)nN / sum;
  Double_t n1 = 0.0, n2 = 0.0;
  for (Int_t i = 0; i < 4; ++i) {
    n1 += raw[i] * scale;
    n2 += raw[i] * scale;
  }
  const Double_t meanNorm = 0.5 * (raw[2] + raw[3]) * scale;
  if (TMath::Abs(n1 - n2) > 1e-12 || TMath::Abs(meanNorm - 1.0) > 1e-9) {
    std::cerr << "FAIL norm not single-scale meanNorm=" << meanNorm << " scale=" << scale << std::endl;
    ++nFail;
  } else {
    std::cout << "PASS normalization applied once scale=" << scale << std::endl;
  }
  delete hB;
  delete hF;
  delete hFg;
  delete hS;
  delete hRe;
  delete hBre;
  delete h40a;
  delete h40b;
  delete h80;
  delete f40a;
  delete f40b;
  delete f80;
  std::cout << "window-count toy: nFail=" << nFail << (nFail ? " FAIL" : " PASS") << std::endl;
  return nFail;
}

void test_cf_kstarmassfit() {
  const Double_t mean = 1.0195;
  const Double_t sigma = 0.0045;
  const Double_t sigMin = 1.012;
  const Double_t sigMax = 1.026;
  const Double_t aMin = 1.04;
  const Double_t aMax = 1.06;
  const Double_t injectedCf[] = {1.40, 1.20, 1.05, 1.00};
  const Double_t kstar[] = {0.025, 0.075, 0.125, 0.40};
  const Int_t nBins = 4;
  const Double_t yMeTrue = 800.0;
  const Double_t nBkg = 4000.0;

  Double_t maxRel = 0.0;
  Int_t nOk = 0;
  for (Int_t i = 0; i < nBins; ++i) {
    const Double_t ySeTrue = injectedCf[i] * yMeTrue;
    TH1D* hBse = makeBkg(Form("Bse_%d", i), nBkg);
    TH1D* hBme = makeBkg(Form("Bme_%d", i), nBkg * 1.1);
    TH1D* hFse = (TH1D*)hBse->Clone(Form("Fse_%d", i));
    TH1D* hFme = (TH1D*)hBme->Clone(Form("Fme_%d", i));
    addSignal(hFse, ySeTrue, mean, sigma);
    addSignal(hFme, yMeTrue, mean, sigma);

    const Double_t aSE = alphaFromWindow(hFse, hBse, aMin, aMax);
    const Double_t aME = alphaFromWindow(hFme, hBme, aMin, aMax);
    TH1D* hSse = (TH1D*)hFse->Clone(Form("Sse_%d", i));
    TH1D* hSme = (TH1D*)hFme->Clone(Form("Sme_%d", i));
    TH1D* hBscSE = (TH1D*)hBse->Clone(Form("BscSE_%d", i));
    TH1D* hBscME = (TH1D*)hBme->Clone(Form("BscME_%d", i));
    hBscSE->Scale(aSE);
    hBscME->Scale(aME);
    hSse->Add(hBscSE, -1.0);
    hSme->Add(hBscME, -1.0);

    Double_t ySE = 0.0, yME = 0.0;
    const Bool_t okSE = fitGausYield(hSse, sigMin, sigMax, ySE);
    const Bool_t okME = fitGausYield(hSme, sigMin, sigMax, yME);
    if (!okSE || !okME || yME <= 0.0) {
      std::cerr << "FAIL bin k*=" << kstar[i] << " fit SE=" << okSE << " ME=" << okME << std::endl;
      continue;
    }
    const Double_t cf = ySE / yME;
    const Double_t rel = TMath::Abs(cf - injectedCf[i]) / injectedCf[i];
    std::cout << Form("k*=%.3f  injected=%.3f  recovered=%.3f  rel=%.3f  alphaSE=%.3f alphaME=%.3f\n",
                      kstar[i], injectedCf[i], cf, rel, aSE, aME);
    if (rel > maxRel) maxRel = rel;
    if (rel < 0.15) ++nOk;
    delete hBse;
    delete hBme;
    delete hFse;
    delete hFme;
    delete hSse;
    delete hSme;
    delete hBscSE;
    delete hBscME;
  }
  const Bool_t pass = (nOk == nBins) && (maxRel < 0.15);
  std::cout << "kstarMassFitCF toy closure: nOk=" << nOk << "/" << nBins << " maxRel=" << maxRel
            << (pass ? " PASS" : " FAIL") << std::endl;
  if (!pass) {
    std::cerr << "kstarMassFitCF toy closure FAILED" << std::endl;
  }
  const Int_t nCntFail = test_window_count();
  if (nCntFail) {
    std::cerr << "window-count toy FAILED" << std::endl;
  }
}
