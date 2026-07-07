
#include "Offline/Validation/inc/TValHistE.hh"
#include "TLatex.h"
#include "TMath.h"
#include "TPad.h"
#include "TStyle.h"

// ClassImp(TValHistE)

//_____________________________________________________________________________
void TValHistE::Clear(Option_t* Opt) {
  fEff1 = NULL;
  fEff2 = NULL;
  fSum1 = 0.0;
  fSum2 = 0.0;
  ClearB();
}

//_____________________________________________________________________________
Int_t TValHistE::Analyze(Option_t* Opt) {
  fKsProb = 0.0;
  fFrProb = 0.0;
  fDiff = true;
  fStatus = fCantCompare;

  if (fEff1 == NULL || fEff2 == NULL) {
    return fStatus;
  }

  TString name1 = fEff1->ClassName();
  if (name1 != "TEfficiency") {
    return fStatus;
  }
  TString name2 = fEff2->ClassName();
  if (name2 != "TEfficiency") {
    return fStatus;
  }

  if (fEff1->GetTotalHistogram()->GetNbinsX() <= 0 ||
      fEff1->GetTotalHistogram()->GetNbinsX() !=
          fEff2->GetTotalHistogram()->GetNbinsX()) {
    return fStatus;
  }

  uint llim = (fPar.GetUnder() == 0 ? 0 : 1);
  uint ulim =
      fEff1->GetTotalHistogram()->GetNbinsX() + (fPar.GetOver() == 0 ? 1 : 0);

  // do sums
  fSum1 = 0.0;
  fSum2 = 0.0;
  fDiff = false;

  for (uint ii = llim; ii <= ulim; ii++) {
    // if not identical, set the flag
    if (fEff1->GetEfficiency(ii) != fEff2->GetEfficiency(ii)) fDiff = true;
    fSum1 += fEff1->GetTotalHistogram()->GetBinContent(ii);
    fSum2 += fEff2->GetTotalHistogram()->GetBinContent(ii);
  }

  // if both zero, they agree
  if (fSum1 == 0 && fSum2 == 0) {
    fEmpty = true;
    fKsProb = 1.0;
    fFrProb = 1.0;
    fStatus = fPerfect;
    return fStatus;
  }

  // no sensible meaning to fractional difference
  fFrProb = 1.0;

  if (fSum1 == 0 || fSum2 == 0) {
    fEmpty = true;
    // prob of observing zero event when expecting N
    fKsProb = TMath::Exp(-TMath::Max(fSum1, fSum2));
    fFrProb = 0.0;
  } else {
    // KsProb is filled from a chi2 comparison
    int ndof = 0;
    double chi2 = 0.0;
    double n1, n2, c1, c2, e1, e2;
    for (uint ii = llim; ii <= ulim; ii++) {
      n1 = fEff1->GetTotalHistogram()->GetBinContent(ii);
      n2 = fEff2->GetTotalHistogram()->GetBinContent(ii);
      c1 = fEff1->GetEfficiency(ii);
      c2 = fEff2->GetEfficiency(ii);
      e1 =
          (fEff1->GetEfficiencyErrorLow(ii) + fEff1->GetEfficiencyErrorUp(ii)) /
          2.0;
      e2 =
          (fEff2->GetEfficiencyErrorLow(ii) + fEff2->GetEfficiencyErrorUp(ii)) /
          2.0;
      if (n1 > 0 && n2 > 0) {
        // both have entries, add to chi2
        chi2 += pow(c1 - c2, 2) / (pow(e1, 2) + pow(e2, 2));
        ndof++;
      } else if (n1 > 0 || n2 > 0) {
        // one has entries, the other doesn't, we need some penalty
        // the idea is to add the Poisson probability of observing
        // zero when expecting n to the likelihood, represented by
        // adding n to the chi2
        chi2 += (n1 > 0 ? n1 : n2);
        ndof++;
      }
      // if both effs have no entries, ignore this bin

    }  // loop over bins

    if (ndof > 0) {
      fKsProb = TMath::Prob(chi2, ndof);
    } else {
      // shouldn't come here since we returned above if no entries
      // in either plot.  If there were entries then ndof should be >0
      fKsProb = 0.0;
    }
  }  // end if a sum is zero

  fStatus = fFail;
  if (fKsProb > fPar.GetLoose()) fStatus = fLoose;
  if (fKsProb > fPar.GetTight()) fStatus = fTight;
  if (!fDiff) fStatus = fPerfect;

  return fStatus;
}

//_____________________________________________________________________________
void TValHistE::Summary(Option_t* Opt) {
  printf("%8.5f %8.5f %2d %10g %10g %s/%s \"%s\"\n", fKsProb, fFrProb, fStatus,
         fSum1, fSum2, GetTag().Data(), GetName(), GetTitle());
}

//_____________________________________________________________________________
void TValHistE::Draw(Option_t* Opt) {
  TString opt1 = Opt;
  opt1.ToLower();
  // TEfficiency doesn't handle log scale option, so ignore it
  bool qLog = false;

  if (fEff1->GetDimension() != 1) {
    printf("Draw for 2D efficiency is not implemented\n");
    return;
  }

  // int color1 = kRed+1;
  // int color2 = kGreen+2;
  // int color1 = kGray+1;
  // int color2 = kMagenta+2;
  // int color1 = kAzure-4;
  // int color2 = kRed;
  // int color1 = kGreen-9;
  int color1 = kBlack;
  int color2 = kRed + 1;

  double x;
  double xmax = 0.;
  double xmin = 1.0e20;
  double xmnz = 1.0e20;  // min not counting 0, for log scales

  int nx = fEff1->GetTotalHistogram()->GetNbinsX();
  for (int i = 1; i <= nx; i++) {
    x = fEff1->GetEfficiency(i) + fEff2->GetEfficiencyErrorUp(i);
    if (x > xmax) xmax = x;
    x = fEff2->GetEfficiency(i) + fEff2->GetEfficiencyErrorUp(i);
    if (x > xmax) xmax = x;
    x = fEff1->GetEfficiency(i) - fEff1->GetEfficiencyErrorLow(i);
    if (x < xmin) xmin = x;
    if (x < xmnz && x > 0.0) xmnz = x;
    x = fEff2->GetEfficiency(i) - fEff2->GetEfficiencyErrorLow(i);
    if (x < xmin) xmin = x;
    if (x < xmnz && x > 0.0) xmnz = x;
  }

  uint llim = (fPar.GetUnder() == 0 ? 0 : 1);
  uint ulim =
      fEff1->GetTotalHistogram()->GetNbinsX() + (fPar.GetOver() == 0 ? 1 : 0);
  double ymean1 = 0.0;
  double ymean2 = 0.0;
  int n1 = 0, n2 = 0;
  for (uint ii = llim; ii <= ulim; ii++) {
    if (fEff1->GetTotalHistogram()->GetBinContent(ii) > 0) {
      ymean1 += fEff1->GetEfficiency(ii);
      n1++;
    }
    if (fEff2->GetTotalHistogram()->GetBinContent(ii) > 0) {
      ymean2 += fEff2->GetEfficiency(ii);
      n2++;
    }
  }
  ymean1 = (n1 > 0 ? ymean1 / n1 : 0.0);
  ymean2 = (n2 > 0 ? ymean2 / n2 : 0.0);

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0);

  if (qLog) {
    gPad->SetLogy(1);
  } else {
    gPad->SetLogy(0);
  }

  fEff1->SetMarkerSize(fFontScale * 0.7);
  fEff2->SetMarkerStyle(21);
  fEff1->SetMarkerColor(color1);
  fEff1->SetLineColor(color1);
  fEff1->Draw("");

  fEff2->SetMarkerSize(fFontScale * 0.6);
  fEff2->SetMarkerStyle(20);
  fEff2->SetMarkerColor(color2);
  fEff2->SetLineColor(color2);
  fEff2->Draw("SAME");

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
  r = (fSum1 != 0 ? fSum2 / fSum1 : 0.0);
  snprintf(tstring, sizeof(tstring), "N %10d %10d %6f   U %10g %10g", int(fSum1), int(fSum2), r,
          fEff1->GetTotalHistogram()->GetBinContent(0),
          fEff2->GetTotalHistogram()->GetBinContent(0));
  tr1->SetNDC();
  tr1->SetText(0.13, 0.87, tstring);
  tr1->SetTextSize(fFontScale * 0.030);
  tr1->SetTextFont(102);
  tr1->Draw("SAME");

  TText* tr2 = new TText();
  r = (ymean1 != 0.0 ? ymean2 / ymean1 : 0.0);
  int iover = fEff1->GetTotalHistogram()->GetNbinsX() + 1;
  snprintf(tstring, sizeof(tstring), "Y %10g %10g %6f   O %10g %10g", ymean1, ymean2, r,
          fEff1->GetTotalHistogram()->GetBinContent(iover),
          fEff2->GetTotalHistogram()->GetBinContent(iover));
  tr2->SetNDC();
  tr2->SetText(0.13, 0.84, tstring);
  tr2->SetTextSize(fFontScale * 0.030);
  tr2->SetTextFont(102);
  tr2->Draw();

  /*

  TText* tr3 = new TText();
  r = (fEff1->GetEntries()>0.0?
       fEff2->GetEntries()/fEff1->GetEntries() : 1.0);
  snprintf(tstring, sizeof(tstring), "N %10g %10g %6f",
          fEff1->GetEntries(),fEff2->GetEntries(),r);
  tr3->SetNDC();
  tr3->SetText(0.13,0.81,tstring);
  tr3->SetTextSize(fFontScale*0.030);
  tr3->SetTextFont(102);
  tr3->Draw();
  */

  gPad->Update();
}

//_____________________________________________________________________________
void TValHistE::Dump() const {
  TObject::Dump();
  printf("bin        x         Hist1        Hist2       Diff\n");
  if (fEff1 && fEff2) {
    for (int i = 0; i <= fEff1->GetTotalHistogram()->GetNbinsX() + 1; i++) {
      printf("%3d %10.5f %10.0f %10.0f %10.0f\n", i,
             fEff1->GetTotalHistogram()->GetBinCenter(i),
             fEff1->GetEfficiency(i), fEff2->GetEfficiency(i),
             fEff1->GetEfficiency(i) - fEff2->GetEfficiency(i));
    }
  }
}
