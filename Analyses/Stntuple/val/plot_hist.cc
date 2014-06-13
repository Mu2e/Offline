///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "TH1.h"
#include "TText.h"
#include "TInterpreter.h"
#include "TSystem.h"

#include "Stntuple/val/stntuple_val_functions.hh"

void plot_hist(const char* fn, 
	       const char* Module, 
	       const char* Hist, 
	       const char* Opt,
	       const int   MarkerStyle,
	       const int   MarkerSize ,
	       const int   PrintFileName) 
{

  float xmin, xmax, ymax, ymin, x0, y0;
//-----------------------------------------------------------------------------
//  get histogram and make sure it is normalized properly
//-----------------------------------------------------------------------------
  TH1F* h1  = get_hist(fn,Module,Hist);
  //  h1->SetNormFactor(h1->GetEntries());
  h1->SetFillStyle(0);

//    h1->SetStatFontSize(0.04);
//    h1->SetStatFont(42);

  h1->GetXaxis()->SetTitleSize(0.05);
  h1->GetXaxis()->SetTitleFont(42);
  h1->GetXaxis()->SetTitleOffset(0.7);

  if (MarkerStyle > 0) h1->SetMarkerStyle(MarkerStyle);
  if (MarkerSize  > 0) h1->SetMarkerSize (MarkerSize );

  h1->Draw(Opt);

  
  if (PrintFileName) {
    xmin = h1->GetXaxis()->GetXmin();
    xmax = h1->GetXaxis()->GetXmax();

    if (h1->InheritsFrom("TH2") == 0) {
      ymin = 0;
      ymax = h1->GetMaximum();
    }
    else {
      ymin = h1->GetYaxis()->GetXmin();
      ymax = h1->GetYaxis()->GetXmax();
    }

    x0 = xmin+(xmax-xmin)/10;
    y0 = ymax-(ymax-ymin)/20;

    TText* text = new TText(x0,y0,fn);
    text->SetTextFont(42);
    text->SetTextSize(0.03);
    text->Draw();
  }

}
