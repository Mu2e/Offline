///////////////////////////////////////////////////////////////////////////////
// Nov 30 2001 P.Murat
// plot turn-on curve of L1 MET25 trigger
// bin 5 GeV 0-100 GeV
//
// input: met25_ana_hist.root
//
// 
// met  20     25     30    35     40   45 
// eff  0.18   0.44   0.73  0.92   1    1
///////////////////////////////////////////////////////////////////////////////

#if !defined (__CINT__) || defined (__MAKECINT__)
#include <cmath>

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TFile.h"
#include "TFolder.h"
#include "Stntuple/val/stntuple_val_functions.hh"

#endif

TH1F* create_eff_hist(const char* FName,
		      const char* Module,
		      const char* DenomHist, 
		      const char* NumerHist,
		      const char* EffHist,
		      const char* EffTitle,
		      Int_t       MarkerStyle, 
		      Float_t     MarkerSize) 
{

//   const char* fname = 
//     "/usr/people/murat/www/tau/wtau/monojet/hist/met25_ana_hist.root";

//-----------------------------------------------------------------------------
//  make sure all the needed scripts are loaded
//-----------------------------------------------------------------------------
  // if (! gInterpreter->IsLoaded("plot/get_hist.C")) {
  //   gInterpreter->LoadMacro("plot/get_hist.C");
  // }

  TH1F* denom = gh1(FName,Module,DenomHist);
  TH1F* numer = gh1(FName,Module,NumerHist);

  TH1F* eff   = (TH1F*) numer->Clone(EffHist);

  if (strcmp(EffTitle,"") == 0) eff->SetTitle("Efficiency");
  else                          eff->SetTitle(EffTitle);

  int nbins = eff->GetNbinsX();

  double err;

  for (int i=1; i<=nbins; i++) {
    double n1 = numer->GetBinContent(i);
    double n2 = denom->GetBinContent(i);

    if (n2 > 0) {
      eff->SetBinContent(i,n1/n2);
      if (n1/n2 < 0.5) {
	err = sqrt(n1+1.e-6)/n2;
      }
      else {
	err = sqrt(n2-n1+1.e-6)/n2;
      }

      eff->SetBinError(i,err);
    }
    else {
      eff->SetBinContent(i,0);
      eff->SetBinError(i,0);
    }
  }
  eff->SetMarkerStyle(MarkerStyle);
  eff->SetMarkerSize (MarkerSize);
  eff->SetMinimum(0);
  eff->SetMaximum(1.5);

  return eff;
}
