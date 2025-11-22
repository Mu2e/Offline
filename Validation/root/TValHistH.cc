
#include "Offline/Validation/inc/TValHistH.hh"
#include "TLatex.h"
#include "TMath.h"
#include "TPad.h"
#include "TStyle.h"

// ClassImp(TValHistH)

//_____________________________________________________________________________
void TValHistH::Clear(Option_t* Opt) {
  fHist1 = NULL;
  fHist2 = NULL;
  fSum1 = 0.0;
  fSum2 = 0.0;
  fNorm2 = 0.0;
  ClearB();
}

//_____________________________________________________________________________
Int_t TValHistH::Analyze(Option_t* Opt) {
  fKsProb = 0.0;
  fFrProb = 0.0;
  fDiff = true;
  fStatus = fCantCompare;

  if (fHist1 == NULL || fHist2 == NULL) {
    return fStatus;
  }

  TString name1 = fHist1->ClassName();
  if (name1 != "TH1F" && name1 != "TH1D") {
    return fStatus;
  }
  TString name2 = fHist2->ClassName();
  if (name2 != "TH1F" && name2 != "TH1D") {
    return fStatus;
  }

  if (fHist1->GetNbinsX() <= 0 || fHist1->GetNbinsX() != fHist2->GetNbinsX()) {
    return fStatus;
  }

  uint llim = (fPar.GetUnder() == 0 ? 0 : 1);
  uint ulim = fHist1->GetNbinsX() + (fPar.GetOver() == 0 ? 1 : 0);

  // do sums
  fSum1 = 0.0;
  fSum2 = 0.0;
  fDiff = false;

  for (uint ii = llim; ii <= ulim; ii++) {
    // if not identical, set the flag
    if (fHist1->GetBinContent(ii) != fHist2->GetBinContent(ii)) fDiff = true;
    fSum1 += fHist1->GetBinContent(ii);
    fSum2 += fHist2->GetBinContent(ii);
  }

  // find normalization for hist 2
  if (fSum2 > 0) {
    double t2 = fHist2->Integral();  // count ignoring under/over
    // root SetNormFactor applies to non-under/overflow sum
    if (fPar.GetMode() == 1) {
      fNorm2 = (fSum1 / fSum2) * t2;
    } else if (fPar.GetMode() == 2 && fPar.GetScale2() > 0.0) {
      // scale set by user
      fNorm2 = (fPar.GetScale1() / fPar.GetScale2()) * t2;
    }
  }
  if (fNorm2 > 0.0) fHist2->SetNormFactor(fNorm2);

  // if both zero, they agree
  if (fSum1 == 0 && fSum2 == 0) {
    fEmpty = true;
    fKsProb = 1.0;
    fFrProb = 1.0;
    fStatus = fPerfect;
    return fStatus;
  }

  // find fractional difference
  if (fSum1 == 0 || fSum2 == 0) {
    fEmpty = true;
    fFrProb = 0.0;
  } else {
    Double_t maxDiff = 0.0;
    Double_t s1 = 0.0, s2 = 0.0;
    for (uint ii = llim; ii <= ulim; ii++) {
      s1 += fPar.GetScale1() * fHist1->GetBinContent(ii);
      s2 += fPar.GetScale2() * fHist2->GetBinContent(ii);
      if (TMath::Abs(s1 - s2) > maxDiff) maxDiff = TMath::Abs(s1 - s2);
    }
    fFrProb = 1.0 - maxDiff / s1;  // 1 is the standard
    if (fFrProb < 0.0) fFrProb = 0.0;
  }

  if (fSum1 == 0 || fSum2 == 0) {
    // prob of observing zero event when expecting N
    fKsProb = TMath::Exp(-TMath::Max(fSum1, fSum2));
  } else {
    TString qual;
    if (fPar.GetUnder() == 0) qual.Append("U");
    if (fPar.GetOver() == 0) qual.Append("O");
    fKsProb = fHist1->KolmogorovTest(fHist2, qual.Data());
  }

  fStatus = fFail;
  if (fPar.GetIndependent() == 0) {
    if (fFrProb > fPar.GetLoose() || fKsProb > fPar.GetLoose())
      fStatus = fLoose;
    if (fFrProb > fPar.GetTight() || fKsProb > fPar.GetTight())
      fStatus = fTight;
  } else {
    if (fKsProb > fPar.GetLoose()) fStatus = fLoose;
    if (fKsProb > fPar.GetTight()) fStatus = fTight;
  }
  if (!fDiff) fStatus = fPerfect;

  return fStatus;
}

//_____________________________________________________________________________
void TValHistH::Summary(Option_t* Opt) {
  printf("%8.5f %8.5f %2d %10g %10g %s/%s \"%s\"\n", fKsProb, fFrProb, fStatus,
         fSum1, fSum2, GetTag().Data(), GetName(), GetTitle());
}

//_____________________________________________________________________________
void TValHistH::Draw(Option_t* Opt) {
  TString opt1 = Opt;
  opt1.ToLower();
  bool qLog = (opt1.Index("log") >= 0);

  if (fHist1->InheritsFrom("TH2")) {
    printf("Draw for 2D is not implemented\n");
    return;
  }

  // int color1 = kRed+1;
  // int color2 = kGreen+2;
  // int color1 = kGray+1;
  // int color2 = kMagenta+2;
  // int color1 = kAzure-4;
  // int color2 = kRed;
  // int color1 = kGreen-9;
  int color1 = kGray;
  int color2 = kRed + 1;

  double x;
  double xmax = 0.;
  double xmin = 1.0e20;
  double xmnz = 1.0e20;  // min not counting 0, for log scales

  int nx = fHist1->GetNbinsX();
  double tot2 = fHist2->Integral();
  double scale2 = ((tot2 > 0.0 && fNorm2 > 0.0) ? fNorm2 / tot2 : 1.0);
  for (int i = 1; i <= nx; i++) {
    x = fHist1->GetBinContent(i);
    if (x > xmax) xmax = x;
    if (x < xmin) xmin = x;
    if (x < xmnz && x > 0.0) xmnz = x;
    x = scale2 * fHist2->GetBinContent(i);
    if (x > xmax) xmax = x;
    if (x < xmin) xmin = x;
    if (x < xmnz && x > 0.0) xmnz = x;
    x = scale2 * (fHist2->GetBinContent(i) + fHist2->GetBinError(i));
    if (x > xmax) xmax = x;
    x = scale2 * (fHist2->GetBinContent(i) - fHist2->GetBinError(i));
    if (x < xmin) xmin = x;
    if (x < xmnz && x > 0.0) xmnz = x;
  }

  // int optstat = gStyle->GetOptStat();

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  if (qLog) {
    float xh = (xmax > 1.0e-20 ? xmax : 1.0);
    float xl = (xmnz < 1.0e19 ? xmnz : 0.1);
    fHist1->SetMaximum(xh * TMath::Power(xh / xl, 0.14));
    fHist1->SetMinimum(xl / TMath::Power(xh / xl, 0.05));
    // if 2's min is 0, log scale will give an error, so set it up too
    fHist2->SetMinimum(fHist1->GetMinimum());
  } else {
    float xh = (xmax > 1.0e-20 ? xmax : 1.0);
    fHist1->SetMaximum(xh * 1.14);
    fHist1->SetMinimum(0.0);
  }

  fHist1->SetFillColor(color1);
  fHist1->SetLineColor(color1);
  fHist1->Draw("HIST");

  if (qLog) {
    gPad->SetLogy(1);
  } else {
    gPad->SetLogy(0);
  }

  fHist2->SetMarkerSize(fFontScale * 0.6);
  fHist2->SetMarkerStyle(20);
  fHist2->SetMarkerColor(color2);
  fHist2->SetLineColor(kGray + 2);
  fHist2->Draw("e0SAME");

  char tstring[200];
  int color;

  color = kRed;
  if (GetKsProb() > fPar.GetLoose()) color = kOrange;
  if (GetKsProb() > fPar.GetTight()) color = kGreen;
  if (GetStatus() == fPerfect) color = kGreen + 2;
  snprintf(tstring, sizeof(tstring), "KS=%8.6f", GetKsProb());
  TText* t1 = new TText();
  t1->SetNDC();
  t1->SetText(0.15, 0.91, tstring);
  t1->SetTextSize(fFontScale * 0.035);
  t1->SetTextFont(42);
  t1->SetTextColor(color);
  t1->Draw();

  color = kRed;
  if (GetFrProb() > fPar.GetLoose()) color = kOrange;
  if (GetFrProb() > fPar.GetTight()) color = kGreen;
  if (GetStatus() == fPerfect) color = kGreen + 2;
  snprintf(tstring, sizeof(tstring), "FR=%8.6f", GetFrProb());
  TText* t2 = new TText();
  t2->SetNDC();
  t2->SetText(0.32, 0.91, tstring);
  t2->SetTextSize(fFontScale * 0.035);
  t2->SetTextFont(42);
  t2->SetTextColor(color);
  t2->Draw();

  TText* ts = new TText();
  snprintf(tstring, sizeof(tstring), "%s/%s \"%s\"", GetTag().Data(), GetName(), GetTitle());
  ts->SetNDC();
  ts->SetText(0.10, 0.95, tstring);
  ts->SetTextSize(fFontScale * 0.035);
  ts->SetTextFont(42);
  ts->Draw();

  double r;
  TText* tr1 = new TText();
  r = (fHist1->GetMean() != 0.0 ? fHist2->GetMean() / fHist1->GetMean() : 1.0);
  snprintf(tstring, sizeof(tstring), "M %10g %10g %6f   U %10g %10g", fHist1->GetMean(),
          fHist2->GetMean(), r, fHist1->GetBinContent(0),
          fHist2->GetBinContent(0));
  tr1->SetNDC();
  tr1->SetText(0.13, 0.87, tstring);
  tr1->SetTextSize(fFontScale * 0.030);
  tr1->SetTextFont(102);
  tr1->Draw("SAME");

  TText* tr2 = new TText();
  r = (fHist1->GetRMS() > 0.0 ? fHist2->GetRMS() / fHist1->GetRMS() : 1.0);
  int iover = fHist1->GetNbinsX() + 1;
  snprintf(tstring, sizeof(tstring), "R %10g %10g %6f   O %10g %10g", fHist1->GetRMS(),
          fHist2->GetRMS(), r, fHist1->GetBinContent(iover),
          fHist2->GetBinContent(iover));
  tr2->SetNDC();
  tr2->SetText(0.13, 0.84, tstring);
  tr2->SetTextSize(fFontScale * 0.030);
  tr2->SetTextFont(102);
  tr2->Draw();

  TText* tr3 = new TText();
  r = (fHist1->GetEntries() > 0.0 ? fHist2->GetEntries() / fHist1->GetEntries()
                                  : 1.0);
  snprintf(tstring, sizeof(tstring), "N %10g %10g %6f", fHist1->GetEntries(),
          fHist2->GetEntries(), r);
  tr3->SetNDC();
  tr3->SetText(0.13, 0.81, tstring);
  tr3->SetTextSize(fFontScale * 0.030);
  tr3->SetTextFont(102);
  tr3->Draw();

  // Summary();
  gPad->Update();

  /*
  printf(" ------------ %s , ks(prob) : %12.8g\n",GetName(),GetKsProb());
  printf(" xmax %f\n",xmax);
  printf(" h1i, h2i = %f  %f\n",fHist1->Integral(), fHist2->Integral());
  */
  // gStyle->SetOptStat(optstat);
}

//_____________________________________________________________________________
void TValHistH::Dump() const {
  TObject::Dump();
  printf("bin        x         Hist1        Hist2       Diff\n");
  if (fHist1 && fHist2) {
    for (int i = 0; i <= fHist1->GetNbinsX() + 1; i++) {
      printf("%3d %10.5f %10.0f %10.0f %10.0f\n", i, fHist1->GetBinCenter(i),
             fHist1->GetBinContent(i), fHist2->GetBinContent(i),
             fHist1->GetBinContent(i) - fHist2->GetBinContent(i));
    }
  }
}
