///////////////////////////////////////////////////////////////////////////////
// the 1st histogram is drawn with "ep" (data)
// the second one    - as a histogram   (MC)
// assume that we have module_name/Hist/hist_name
///////////////////////////////////////////////////////////////////////////////
#include "TString.h"
#include "TH1.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "Stntuple/val/stntuple_val_functions.hh"

//_____________________________________________________________________________
void overlay_2h(TH1*        h1, 
		TH1*        h2, 
		TString     Opt,
		int         MinXBin,
		int         MaxXBin,
		const float MarkerSize,
		int         Scale     )
{
//-----------------------------------------------------------------------------
//  make sure all the points are displayed
//-----------------------------------------------------------------------------
  double x, scalef;
  double xmax = 0.;

  if (MinXBin == 0) MinXBin = 1;
  if (MaxXBin == 0) MaxXBin = h1->GetNbinsX();

  h1->SetNormFactor(h1->GetEntries());
  h2->SetNormFactor(h2->GetEntries());

  double h1_int = h1->Integral(MinXBin,MaxXBin);
  double h2_int = h2->Integral(MinXBin,MaxXBin);

  printf(" xmax %f\n",xmax);

  if (Scale) {
    scalef = h1_int/h2_int;
    printf("scalefactor = %10.3f\n",scalef);

    int nx = h2->GetNbinsX();
    for (int i=1; i<=nx; i++) {
      x = h1->GetBinContent(i);
      if (x > xmax) xmax = x;
      x = h2->GetBinContent(i)*scalef;
      if (x > xmax) xmax = x;
    }
    h1->SetMaximum(xmax*1.1);
  }

  h1->GetXaxis()->SetRange(MinXBin,MaxXBin+1);

  h1->SetMarkerSize(MarkerSize);
  h1->SetMarkerStyle(20);
  h1->Draw(Opt);

  if (Scale) {
    h1->SetNormFactor(h1_int);
    h2->SetNormFactor(h1_int/h2_int*h2->Integral());
    printf(" h2 norm factor: %f\n",h2->GetNormFactor());
  }

  h2->SetFillStyle(3001);
  h2->SetFillColor(20);
  h2->Draw("same");

  h1->SetMarkerSize(MarkerSize);
  h1->SetMarkerStyle(20);

  TString oopt(Opt);
  oopt += ",same";
  h1->Draw(oopt.Data());


  printf(" xmax %f\n",xmax);
  printf(" h1i, h2i = %f  %f\n",h1->Integral(), h2->Integral());
}


//_____________________________________________________________________________
void overlay_2h_2f(TString     FName1, 
		   TString     FName2, 
		   TString     ModuleName,
		   TString     Hist1Name,
		   TString     Hist2Name,
		   TString     Opt      ,
		   int         MinXBin  ,
		   int         MaxXBin  ,
		   const float MarkerSize,
		   int         Scale     )
{
//-----------------------------------------------------------------------------
//  get histograms
//-----------------------------------------------------------------------------
  TH1F* h1  = get_hist(FName1,ModuleName,Hist1Name);

  TString name = Hist2Name;
  if (name == "") name = Hist1Name;
  TH1F* h2  = get_hist(FName2,ModuleName,name);

  overlay_2h(h1,h2,Opt,MinXBin,MaxXBin,MarkerSize,Scale);
}

//_____________________________________________________________________________
void overlay_2h_1f(TString     FName1, 
		   TString     ModuleName,
		   TString     Hist1Name,
		   TString     Hist2Name,
		   TString     Opt       ,
		   const int   MinXBin   ,
		         int   MaxXBin   ,
		   const float MarkerSize,
		   Int_t       Scale     )
{
  // overlay 2 histograms produced by the same module (from the same file)

//-----------------------------------------------------------------------------
//  get histograms
//-----------------------------------------------------------------------------
  TH1F* h1  = get_hist(FName1,ModuleName,Hist1Name);
  TH1F* h2  = get_hist(FName1,ModuleName,Hist2Name);

  overlay_2h(h1,h2,Opt,MinXBin,MaxXBin,MarkerSize,Scale);
}
