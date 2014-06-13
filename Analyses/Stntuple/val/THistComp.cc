#include "TPad.h"
#include "Stntuple/val/THistComp.hh"
#include "TLatex.h"

ClassImp(THistComp)
ClassImp(TGoodHistComp)
ClassImp(TBadHistComp)

//_____________________________________________________________________________
void THistComp::Draw(Option_t* Opt) {


  if(fHist1->InheritsFrom("TH2")) {
    printf("Draw for 2D is not implemented\n");
    return;
  }

  double x;
  double xmax = 0.;

  double marker_size = 0.5;

  float scale2 = 1.0;
  float sum1 = fHist1->Integral()+fHist1->GetBinContent(0)
    + fHist1->GetBinContent(fHist1->GetNbinsX()+1);
  float sum2 = fHist2->Integral()+fHist2->GetBinContent(0)
    + fHist2->GetBinContent(fHist2->GetNbinsX()+1);

  if(fNorm<0.0 && sum1>0.0 && sum2>0.0) {   // norm 2 to 1's area
    scale2 = sum1/sum2;
  } else if(fNorm>0.0) {   // norm to a fixed ratio
    scale2 = fNorm;
  }
  float norm2 = fHist2->Integral()*scale2;

  int nx = fHist1->GetNbinsX();
  for (int i=1; i<=nx; i++) {
    x = fHist1->GetBinContent(i);
    if (x > xmax) xmax = x;
    x = scale2*fHist2->GetBinContent(i);
    if (x > xmax) xmax = x;
  }

  fHist1->SetMaximum(xmax*1.1);

  fHist1->SetMarkerSize(marker_size);
  fHist1->SetMarkerStyle(20);
  fHist1->Draw(Opt);

  fHist2->SetNormFactor(norm2);

  fHist2->SetFillStyle(3013);
  fHist2->SetFillColor(41);
  fHist2->Draw("same");

  fHist1->SetMarkerSize(marker_size);
  fHist1->SetMarkerStyle(20);

  TString opt(Opt);
  opt += ",same";
  fHist1->Draw(opt.Data());

  TLatex* tt = new TLatex();
  char tstring[200];
  sprintf(tstring,"PROB=%10.8f",GetKsProb());
  tt->SetNDC();
  tt->SetText(0.5,0.90,tstring);
  tt->SetTextSize(0.04);
  tt->Draw();

  gPad->Update();

  printf(" ------------ %s , ks(prob) : %12.8g\n",GetName(),GetKsProb());
  //  printf(" xmax %f\n",xmax);
  printf(" h1i, h2i = %f  %f\n",fHist1->Integral(), fHist2->Integral());
}

//_____________________________________________________________________________
void THistComp::Dump() const {

  TObject::Dump();
  printf("bin        x         Hist1        Hist2       Diff\n");
  if(fHist1 && fHist2) {
    for(int i=0; i<=fHist1->GetNbinsX()+1; i++) {
      printf("%3d %10.5f %10.0f %10.0f %10.0f\n",i,fHist1->GetBinCenter(i),
	     fHist1->GetBinContent(i),
	     fHist2->GetBinContent(i),
	     fHist1->GetBinContent(i)-fHist2->GetBinContent(i));
    }
  }


}
