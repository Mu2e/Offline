#include "TH1F.h"
#include "TF1.h"
#include "TTree.h"
#include "TLegend.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TLine.h"
#include "TArrow.h"
#include "TCut.h"
#include "TBox.h"
#include "TMath.h"
#include "TProfile.h"
#include "TDirectory.h"
#include "Math/Math.h"
#include "THStack.h"
void dio(){

  TCanvas* dioc = new TCanvas("dioc","dio",1200,800);
  dioc->Divide(2,2);

  Double_t dmhi = trueconvmom;
  Double_t dmlow = trueconvmom - diogenrange;
  TH1F* diogen = new TH1F("diogen","True DIO momentum;MeV",nbins,dmlow,dmhi);
  TH1F* evtwt = new TH1F("evtwt","True DIO momentum;MeV",nbins,dmlow,dmhi);
//  evtwt->Sumw2();
  dio->Project("diogen","mcmom");
  dio->Project("evtwt","mcmom","evtwt");
  evtwt->Scale(dioscale);
  evtwt->SetLineColor(kBlue);
  diogen->SetLineColor(kRed);
  evtwt->SetStats(0);
  diogen->SetStats(0);


  Int_t colors[4] = {kRed,kBlue,kGreen,kBlack};
  TH1F* diogenwin[4] = {0,0,0,0};
  TH1F* diodiffwin[4] = {0,0,0,0};
  const char* dopt[4] = {"","same","same","same"};
  const char* cutset[4] = {"Cutset A","Cutset B","Cutset C","Cutset D"};
  TLegend* dgenwinleg = new TLegend(.5,.6,.7,.9);
  for(unsigned icut=0;icut<4;icut++){
    char diogenname[50], diodiffname[50];
    snprintf(diogenname,50,"diogenwin%i",icut);
    diogenwin[icut] = new TH1F(diogenname,"True momentum of DIO in signal box;MeV",100,dmlow,dmhi);
    diogenwin[icut]->SetStats(0);

    snprintf(diodiffname,50,"diodiffwin%i",icut);
    diodiffwin[icut] = new TH1F(diodiffname,"Reco - true momentum of DIO in signal box;MeV",100,-1,2);
    diodiffwin[icut]->SetStats(0);
//   TH1F* diogood[icut] = new TH1F("diogood","True DIO momentum",100,dmlow,dmhi);
//    dio->Project("diogoodwt","mcmom",goodfit);

    TCut quality = ncuts[icut] && t0cuts[icut] && momcuts[icut] && fitcuts[icut];
    TCut final = (reco+pitch+livegate+quality);
    dio->Project(diogenname,"mcent.mom","evtwt"*(final+momwin));
    diogenwin[icut]->SetFillColor(colors[icut]);
    dio->Project(diodiffname,"fit.mom-mcent.mom","evtwt"*(final+momwin));
    diodiffwin[icut]->SetFillColor(colors[icut]);
    dgenwinleg->AddEntry(diogenwin[icut],cutset[icut],"f");
  }


  dioc->cd(1);
  gPad->SetLogy();
// dead-reconing on spectrum, accounting for bins
  double diofscale = ndecay*(dmhi-dmlow)/nbins;
  diocz_f->SetParameter(0,diofscale);
  evtwt->Draw();
  diocz_f->Draw("same");
  diogen->Draw("same");
  TLegend* dioleg = new TLegend(.2,.4,.6,.6);
  dioleg->AddEntry(diogen,"Generated","l");
  dioleg->AddEntry(evtwt,"Weighted","l");
  dioleg->AddEntry(diocz_f,"Czarnecki etal","l");
  dioleg->Draw();

  dioc->cd(2);
  for(unsigned icut=0;icut<4;icut++){
    diospec[icut]->SetFillColor(colors[icut]);
    if(icut==0)
      diospec[icut]->Draw("Hist");
    else
      diospec[icut]->Draw("Histsame");
  }
  dgenwinleg->Draw();
  TLine* momlowl = new TLine(momlow,0.0,momlow,1.5*diospec[0]->GetBinContent(diospec[0]->GetMaximumBin()));
  momlowl->SetLineColor(kBlack);
  momlowl->SetLineStyle(2);
  momlowl->SetLineWidth(2);
  momlowl->Draw();

  TLine* momhighl = new TLine(momhigh,0.0,momhigh,1.5*diospec[0]->GetBinContent(diospec[0]->GetMaximumBin()));
  momhighl->SetLineColor(kBlack);
  momhighl->SetLineStyle(2);
  momhighl->SetLineWidth(2);
  momhighl->Draw();

  dioc->cd(3);
  gPad->SetLogy();
  for(unsigned icut=0;icut<4;icut++){
   diogenwin[icut]->Draw(dopt[icut]);
  }
  dgenwinleg->Draw();

  dioc->cd(4);
  gPad->SetLogy();
  for(unsigned icut=0;icut<4;icut++){
   diodiffwin[icut]->Draw(dopt[icut]);
  }

  dioc->SaveAs((std::string("diocan")+ssuf).c_str());

  TCanvas* diores = new TCanvas("diores","DIO result",800,600);
  gPad->SetLogy();
  diodiffwin[mu2ecut]->Draw();
  double split = 0.35;
  TLine* td = new TLine(split,0.0,split,diodiffwin[mu2ecut]->GetMaximum());
  td->SetLineColor(kBlack);
  td->SetLineStyle(2);
  td->SetLineWidth(2);
  td->Draw();
 
  int istart = diodiffwin[mu2ecut]->FindFixBin(split);
  double core = diodiffwin[mu2ecut]->Integral(0,istart);
  double tail = diodiffwin[mu2ecut]->Integral(istart+1,100);
  double total = core+tail;
  core /= total;
  tail /= total;
  core *= 100;
  tail *= 100;
  cout <<"core = " << core << " tail = " << tail << endl;
  diores->SaveAs("diores.png");
  char ccore[30], ctail[30];
  snprintf(ccore,30,"%4.1f%%",core);
  snprintf(ctail,30,"%4.1f%%",tail);
  TText* tcore = new TText(0.2,0.4,ccore);
  TText* ttail = new TText(0.6,0.4,ctail);
  tcore->SetNDC();
  ttail->SetNDC();
  tcore->Draw();
  ttail->Draw();
}
