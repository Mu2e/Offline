
#include "Offline/Validation/inc/TValHistP.hh"
#include "TLatex.h"
#include "TMath.h"
#include "TPad.h"
#include "TStyle.h"

// ClassImp(TValHistP)

//_____________________________________________________________________________
void TValHistP::Clear(Option_t* Opt) {
  fProf1 = NULL;
  fProf2 = NULL;
  fSum1 = 0.0;
  fSum2 = 0.0;
  ClearB();
}

//_____________________________________________________________________________
Int_t TValHistP::Analyze(Option_t* Opt) {
  fKsProb = 0.0;
  fFrProb = 0.0;
  fDiff = true;
  fStatus = fCantCompare;

  if (fProf1 == NULL || fProf2 == NULL) {
    return fStatus;
  }

  TString name1 = fProf1->ClassName();
  if (name1 != "TProfile") {
    return fStatus;
  }
  TString name2 = fProf2->ClassName();
  if (name2 != "TProfile") {
    return fStatus;
  }

  if (fProf1->GetNbinsX() <= 0 || fProf1->GetNbinsX() != fProf2->GetNbinsX()) {
    return fStatus;
  }

  uint llim = (fPar.GetUnder() == 0 ? 0 : 1);
  uint ulim = fProf1->GetNbinsX() + (fPar.GetOver() == 0 ? 1 : 0);

  // do sums
  fSum1 = 0.0;
  fSum2 = 0.0;
  fDiff = false;
  for (uint ii = llim; ii <= ulim; ii++) {
    // if not identical, set the flag
    if (fProf1->GetBinContent(ii) != fProf2->GetBinContent(ii)) fDiff = true;
    fSum1 += fProf1->GetBinEntries(ii);
    fSum2 += fProf2->GetBinEntries(ii);
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
      n1 = fProf1->GetBinEntries(ii);
      n2 = fProf2->GetBinEntries(ii);
      c1 = fProf1->GetBinContent(ii);
      c2 = fProf2->GetBinContent(ii);
      e1 = fProf1->GetBinError(ii);
      e2 = fProf2->GetBinError(ii);
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
      // if both profs have no entries, ignore this bin

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
void TValHistP::Summary(Option_t* Opt) {
  printf("%8.5f %8.5f %2d %10g %10g %s/%s \"%s\"\n", fKsProb, fFrProb, fStatus,
         fSum1, fSum2, GetTag().Data(), GetName(), GetTitle());
}

//_____________________________________________________________________________
void TValHistP::Draw(Option_t* Opt) {
  TString opt1 = Opt;
  opt1.ToLower();
  bool qLog = (opt1.Index("log") >= 0);

  if (fProf1->GetDimension() != 1) {
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
  int color1 = kBlack;
  int color2 = kRed + 1;

  double x;
  double xmax = 0.;
  double xmin = 1.0e20;
  double xmnz = 1.0e20;  // min not counting 0, for log scales

  int nx = fProf1->GetNbinsX();
  for (int i = 1; i <= nx; i++) {
    x = fProf1->GetBinContent(i) + fProf2->GetBinError(i);
    if (x > xmax) xmax = x;
    x = fProf2->GetBinContent(i) + fProf2->GetBinError(i);
    if (x > xmax) xmax = x;
    x = fProf1->GetBinContent(i) - fProf1->GetBinError(i);
    if (x < xmin) xmin = x;
    if (x < xmnz && x > 0.0) xmnz = x;
    x = fProf2->GetBinContent(i) - fProf2->GetBinError(i);
    if (x < xmin) xmin = x;
    if (x < xmnz && x > 0.0) xmnz = x;
  }

  uint llim = (fPar.GetUnder() == 0 ? 0 : 1);
  uint ulim = fProf1->GetNbinsX() + (fPar.GetOver() == 0 ? 1 : 0);
  double ymean1 = 0.0;
  double ymean2 = 0.0;
  int n1 = 0, n2 = 0;
  for (uint ii = llim; ii <= ulim; ii++) {
    if (fProf1->GetBinEntries(ii) > 0) {
      ymean1 += fProf1->GetBinContent(ii);
      n1++;
    }
    if (fProf2->GetBinEntries(ii) > 0) {
      ymean2 += fProf2->GetBinContent(ii);
      n2++;
    }
  }
  ymean1 = (n1 > 0 ? ymean1 / n1 : 0.0);
  ymean2 = (n2 > 0 ? ymean2 / n2 : 0.0);

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0);

  if (qLog) {
    float xh = (xmax > 1.0e-20 ? xmax : 1.0);
    float xl = (xmnz < 1.0e19 ? xmnz : 0.1);
    fProf1->SetMaximum(xh * TMath::Power(xh / xl, 0.14));
    fProf1->SetMinimum(xl / TMath::Power(xh / xl, 0.05));
  } else {
    fProf1->SetMaximum(xmax + (xmax - xmin) * 0.14);
    fProf1->SetMinimum(xmin - (xmax - xmin) * 0.05);
  }

  if (qLog) {
    gPad->SetLogy(1);
  } else {
    gPad->SetLogy(0);
  }

  fProf1->SetMarkerStyle(21);
  fProf1->SetMarkerSize(fFontScale * 0.7);
  fProf1->SetMarkerColor(color1);
  fProf1->SetLineColor(color1);
  fProf1->Draw("e1");

  fProf2->SetMarkerStyle(20);
  fProf2->SetMarkerSize(fFontScale * 0.6);
  fProf2->SetMarkerColor(color2);
  fProf2->SetLineColor(color2);
  fProf2->Draw("e1SAME");

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
          fProf1->GetBinEntries(0), fProf2->GetBinEntries(0));
  tr1->SetNDC();
  tr1->SetText(0.13, 0.87, tstring);
  tr1->SetTextSize(fFontScale * 0.030);
  tr1->SetTextFont(102);
  tr1->Draw("SAME");

  TText* tr2 = new TText();
  r = (ymean1 != 0.0 ? ymean2 / ymean1 : 0.0);
  int iover = fProf1->GetNbinsX() + 1;
  snprintf(tstring, sizeof(tstring), "Y %10g %10g %6f   O %10g %10g", ymean1, ymean2, r,
          fProf1->GetBinEntries(iover), fProf2->GetBinEntries(iover));
  tr2->SetNDC();
  tr2->SetText(0.13, 0.84, tstring);
  tr2->SetTextSize(fFontScale * 0.030);
  tr2->SetTextFont(102);
  tr2->Draw();

  /*
  TText* tr3 = new TText();
  r = (fProf1->GetEntries()>0.0?
       fProf2->GetEntries()/fProf1->GetEntries() : 1.0);
  snprintf(tstring, sizeof(tstring), "N %10g %10g %6f",
          fProf1->GetEntries(),fProf2->GetEntries(),r);
  tr3->SetNDC();
  tr3->SetText(0.13,0.81,tstring);
  tr3->SetTextSize(fFontScale*0.030);
  tr3->SetTextFont(102);
  tr3->Draw();
  */

  gPad->Update();
}

//_____________________________________________________________________________
void TValHistP::Dump() const {
  TObject::Dump();
  printf("bin        x         Prof1        Prof2       Diff   NSigma\n");
  if (fProf1 && fProf2) {
    for (int i = 0; i <= fProf1->GetNbinsX() + 1; i++) {
      double d = fProf1->GetBinContent(i) - fProf2->GetBinContent(i);
      double s =
          sqrt(pow(fProf1->GetBinError(i), 2) + pow(fProf2->GetBinError(i), 2));
      if (s > 0.0) {
        s = d / s;
      } else {
        s = 0.0;
      }
      printf(
          "%3d %10.5f %10.0f +- %10.0f   %10.0f +- %10.0f   %10.0f  %10.0f\n",
          i, fProf1->GetBinCenter(i), fProf1->GetBinContent(i),
          fProf1->GetBinError(i), fProf2->GetBinContent(i),
          fProf2->GetBinError(i), d, s);
    }
  }
}
