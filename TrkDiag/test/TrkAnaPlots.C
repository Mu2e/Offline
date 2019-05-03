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
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TDirectory.h"
#include "Math/Math.h"
#include "THStack.h"



// The following is from Alexx Perloff, JetMetaAnalysis
double fnc_dscb(double*xx,double*pp) {
  double x   = xx[0];
// gaussian core
  double N   = pp[0];//norm
  double mu  = pp[1];//mean
  double sig = pp[2];//variance
  // transition parameters
  double a1  = pp[3];
  double p1  = pp[4];
  double a2  = pp[5];
  double p2  = pp[6];

  double u   = (x-mu)/sig;
  double A1  = TMath::Power(p1/TMath::Abs(a1),p1)*TMath::Exp(-a1*a1/2);
  double A2  = TMath::Power(p2/TMath::Abs(a2),p2)*TMath::Exp(-a2*a2/2);
  double B1  = p1/TMath::Abs(a1) - TMath::Abs(a1);
  double B2  = p2/TMath::Abs(a2) - TMath::Abs(a2);

  double result(N);
  if      (u<-a1) result *= A1*TMath::Power(B1-u,-p1);
  else if (u<a2)  result *= TMath::Exp(-u*u/2);
  else            result *= A2*TMath::Power(B2+u,-p2);
  return result;
}

void PId(TTree* ta) {

  TH2F* evspep = new TH2F("evspep","Associated Cluster Energy vs Track Momentum;P (MeV/c);E (MeV)",50,0,200,50,0,200);
  TH2F* evspmp = new TH2F("evspmp","Associated Cluster Energy vs Track Momentum;P (MeV/c);E (MeV)",50,0,200,50,0,200);
  TH2F* evspem = new TH2F("evspem","Associated Cluster Energy vs Track Momentum;P (MeV/c);E (MeV)",50,0,200,50,0,200);
  TH2F* evspmm = new TH2F("evspmm","Associated Cluster Energy vs Track Momentum;P (MeV/c);E (MeV)",50,0,200,50,0,200);
  evspep->SetMaximum(50);
  evspmp->SetMaximum(50);
  evspep->SetStats(0);
  evspmp->SetStats(0);
  evspep->SetLineColor(kOrange);
  evspmp->SetLineColor(kCyan);
  evspem->SetMaximum(50);
  evspmm->SetMaximum(50);
  evspem->SetStats(0);
  evspmm->SetStats(0);
  evspem->SetLineColor(kRed);
  evspmm->SetLineColor(kBlue);

  TCanvas* can = new TCanvas("can","can",800,800);

  ta->Project("evspem","dec.eclust:de.mom","de.status>0&&tcnt.ndemc>0&&demc.pdg==11");
  ta->Project("evspmm","dec.eclust:de.mom","de.status>0&&tcnt.ndemc>0&&demc.pdg==13");
  ta->Project("evspep","dec.eclust:de.mom","de.status>0&&tcnt.ndemc>0&&demc.pdg==-11");
  ta->Project("evspmp","dec.eclust:de.mom","de.status>0&&tcnt.ndemc>0&&demc.pdg==-13");

  evspem->Draw();
  evspem->Draw("box");
  evspmm->Draw("boxsame");
  TLegend* leg = new TLegend(0.1,0.6,0.5,0.9);
  leg->AddEntry(evspem,"True e^{-}","L");
  leg->AddEntry(evspmm,"True #mu^{-}","L");
  leg->AddEntry(evspep,"True e^{+}","L");
  leg->AddEntry(evspmp,"True #mu^{+}","L");
  leg->Draw();
}

void MomResp(TTree* ta, double tqcut, double nmu,const char* file="") {
// cuts
  TCut reco("de.status>0");
  char ctext[80];
  snprintf(ctext,80,"detrkqual.trkqual>%f",tqcut);
  TCut goodfit(ctext);
  double tdlow(0.57735027);
  double tdhigh(1.0);
  double t0min(700.0);
  double t0max(1695.0);
  snprintf(ctext,80,"de.td>%5.5f&&de.td<%5.5f",tdlow,tdhigh);
  TCut rpitch = TCut(ctext);
  snprintf(ctext,80,"de.t0>%f&&de.t0<%f",t0min,t0max);
  TCut livegate = TCut(ctext);
  TCut opa = TCut("de.d0<105 && de.d0>-80 && (de.d0+2/de.om)>450 && (de.d0+2/de.om)<680");
  TCut rmomloose("de.mom>100.0");
  TCut physics = rpitch+opa+livegate+rmomloose;

  TF1* dscb = new TF1("dscb",fnc_dscb,-10.0,5,7);
  dscb->SetParName(0,"Norm");
  dscb->SetParName(1,"x0");
  dscb->SetParName(2,"sigma");
  dscb->SetParName(3,"ANeg");
  dscb->SetParName(4,"PNeg");
  dscb->SetParName(5,"APos");
  dscb->SetParName(6,"PPos");

  TCanvas* rcan = new TCanvas("rcan","Momentum Resolution",600,600);
  rcan->Clear();
  gStyle->SetOptFit(111111);
  gStyle->SetOptStat("oumr");
  gPad->SetLogy();
  TH1F* momresp = new TH1F("momresp","momentum response at tracker;MeV/c",251,-10.0,4.0);
  momresp->Sumw2();
  TCut final = reco+goodfit;
//  TCut evtwt = "evtwt.PBIWeight";
//  ta->Project("momresp","de.mom-demcgen.mom",evtwt*final);
  ta->Project("momresp","de.mom-demcgen.mom",final);
  momresp->Scale(1.0/nmu);
  //    ta->Project(mname,"fit.mom-mcent.mom",final);
  double integral = momresp->GetEntries()*momresp->GetBinWidth(1);
  cout << "Integral = " << integral << " mean = " << momresp->GetMean() << " rms = " << momresp->GetRMS() << endl;
  dscb->SetParameters(0.05*integral,-0.6,0.3,0.7,3.0,3.0,3.0);
  dscb->SetNpx(1000);
  dscb->SetParLimits(3,0.0,50.0);
  dscb->SetParLimits(4,1.0,50.0);
  dscb->SetParLimits(5,0.0,50.0);
  dscb->SetParLimits(6,1.0,50.0);

  momresp->SetMinimum(0.5);
  momresp->Fit("dscb","LRQ");
  momresp->Fit("dscb","LRM");
  if(strcmp(file,"")!=0)rcan->SaveAs(file);
}
void MomRes(TTree* ta, double tqcut,double nmu,const char* file="") {
// cuts
  TCut reco("de.status>0");
  char ctext[80];
  snprintf(ctext,80,"detrkqual.trkqual>%f",tqcut);
  TCut goodfit(ctext);
  double tdlow(0.57735027);
  double tdhigh(1.0);
  double t0min(700.0);
  double t0max(1695.0);
  snprintf(ctext,80,"de.td>%5.5f&&de.td<%5.5f",tdlow,tdhigh);
  TCut rpitch = TCut(ctext);
  snprintf(ctext,80,"de.t0>%f&&de.t0<%f",t0min,t0max);
  TCut livegate = TCut(ctext);
  TCut opa = TCut("de.d0<105 && de.d0>-80 && (de.d0+2/de.om)>450 && (de.d0+2/de.om)<680");
  TCut rmomloose("de.mom>100.0");
  TCut physics = rpitch+opa+livegate+rmomloose;

  TF1* dscb = new TF1("dscb",fnc_dscb,-2.0,2.5,7);
  dscb->SetParName(0,"Norm");
  dscb->SetParName(1,"x0");
  dscb->SetParName(2,"sigma");
  dscb->SetParName(3,"ANeg");
  dscb->SetParName(4,"PNeg");
  dscb->SetParName(5,"APos");
  dscb->SetParName(6,"PPos");

  TCanvas* rcan = new TCanvas("rcan","Momentum Resolution",600,600);
  rcan->Clear();
  gStyle->SetOptFit(111111);
  gStyle->SetOptStat("oumr");
  gPad->SetLogy();
  TH1F* momres = new TH1F("momres","momentum resolution at start of tracker;MeV/c",251,-4,4);
  momres->Sumw2();
  TCut final = reco+goodfit+physics;
 // TCut evtwt = "evtwt.PBIWeight";
//  ta->Project("momres","de.mom-demcent.mom",evtwt*final);
  ta->Project("momres","de.mom-sqrt(demcent.momx^2+demcent.momy^2+demcent.momz^2)",final);
  momres->Scale(1.0/nmu);
  //    ta->Project(mname,"fit.mom-mcent.mom",final);
  double integral = momres->GetEntries()*momres->GetBinWidth(1)/nmu;
  cout << "Integral = " << integral << " mean = " << momres->GetMean() << " rms = " << momres->GetRMS() << endl;
  dscb->SetParameters(2*integral,0.0,0.15,1.0,4.5,1.2,10.0);

//  momres->SetMinimum(0.5);
  momres->Fit("dscb","RQ");
  momres->Fit("dscb","RMQ");
  TFitResultPtr fitres = momres->Fit("dscb","SRE");
  cout << "Core Sigma = " << dscb->GetParameter(2) << " +- " << dscb->GetParError(2) << endl;
  cout << "High Side Power = " << dscb->GetParameter(6) << " +- " << dscb->GetParError(6) << endl;

  // count outliers
  int outbin = momres->FindBin(1.1);
  double outint = momres->Integral(outbin,momres->GetNbinsX());
  double totint = momres->Integral();
  double outrat = outint/totint;
  cout <<"Outlier integral = " << outint  << " total integral = " << totint << " Outlier fraction = " << outrat << endl;


  TLine* zero = new TLine(0.0,0.0,0.0,momres->GetBinContent(momres->GetMaximumBin()));
  zero->SetLineStyle(2);
  zero->Draw();

  TPaveText* rtext = new TPaveText(0.1,0.5,0.4,0.9,"NDC");
  rtext->AddText("Reco Cuts");
  char line[40];
  snprintf(line,80,"%4.3f<tan(#lambda)<%4.3f",tdlow,tdhigh);
  rtext->AddText(line);
  snprintf(line,80,"t0>%5.1f nsec",t0min);
  rtext->AddText(line);
  sprintf(line,"%s",goodfit.GetTitle());
  rtext->AddText(line);
  sprintf(line,"%s",rmomloose.GetTitle());
  rtext->AddText(line);
  sprintf(line,"%5.0f Tracks",momres->GetEntries());
  rtext->AddText(line);
  rtext->Draw();
  if(strcmp(file,"")!=0)rcan->SaveAs(file);

}

void Acc(TTree* ta, double tqcut,int ngen,const char* file="") {
  unsigned nbins(8);
  double bmax = nbins-0.5;

  TH1F* acc = new TH1F("acc","CE Acceptance #times Efficiency;;Cummulative a#times#epsilon",nbins,-0.5,bmax);
  TH1F* racc = new TH1F("racc","CE Acceptance #times Efficiency;;Relative a#times#epsilon",nbins,-0.5,bmax);
//  acc->Sumw2();
//  racc->Sumw2();
  unsigned ibin(1);
  acc->GetXaxis()->SetBinLabel(ibin++,"All CE");
  acc->GetXaxis()->SetBinLabel(ibin++,"MC Selection");
  acc->GetXaxis()->SetBinLabel(ibin++,"KF Track fit");
  acc->GetXaxis()->SetBinLabel(ibin++,"Fit Quality");
  acc->GetXaxis()->SetBinLabel(ibin++,"Livegate");
  acc->GetXaxis()->SetBinLabel(ibin++,"Reco pitch");
  acc->GetXaxis()->SetBinLabel(ibin++,"OPA Rejection");
  acc->GetXaxis()->SetBinLabel(ibin++,"Momentum window");


  ibin = 1;
  racc->GetXaxis()->SetBinLabel(ibin++,"All CE");
  racc->GetXaxis()->SetBinLabel(ibin++,"MC Selection");
  racc->GetXaxis()->SetBinLabel(ibin++,"KF Track fit");
  racc->GetXaxis()->SetBinLabel(ibin++,"Fit Quality");
  racc->GetXaxis()->SetBinLabel(ibin++,"Livegate");
  racc->GetXaxis()->SetBinLabel(ibin++,"Reco pitch");
  racc->GetXaxis()->SetBinLabel(ibin++,"OPA Rejection");
  racc->GetXaxis()->SetBinLabel(ibin++,"Momentum window");

  ibin = 0;
  const char* binnames[11] ={"0.0","1.0","2.0","3.0","4.0","5.0","6.0","7.0","8.0","9.0","10.0"};

  TCut mcsel = "demc.ndigigood>15&&sqrt(demcent.momx^2+demcent.momy^2+demcent.momz^2)>90.0";
  // &&demcent.td>0.55&&demcent.td<1.05";
  //&&fmod(demcent.t0,1695.0)>500.0";
  TCut reco = "de.status>0";
  char ctext[80];
  snprintf(ctext,80,"detrkqual.trkqual>%f",tqcut);
  TCut goodfit(ctext);
  TCut livegate = "de.t0>700.0&&de.t0<1695";
  TCut rpitch = "de.td>0.57735027&&de.td<1.0";
  TCut opa = "de.d0<105 && de.d0>-80 && (de.d0+2/de.om)>450 && (de.d0+2/de.om)<680";
  TCut rmom = "de.mom>103.85";
//  TCut evtwt = "evtwt.PBIWeight";
  TCut evtwt = "";
  ta->Project("acc",binnames[ibin++],evtwt);
  ta->Project("+acc",binnames[ibin++],evtwt*mcsel);
  ta->Project("+acc",binnames[ibin++],evtwt*(mcsel+reco));
  ta->Project("+acc",binnames[ibin++],evtwt*(mcsel+reco+goodfit));
  ta->Project("+acc",binnames[ibin++],evtwt*(mcsel+reco+goodfit+livegate));
  ta->Project("+acc",binnames[ibin++],evtwt*(mcsel+reco+goodfit+livegate+rpitch));
  ta->Project("+acc",binnames[ibin++],evtwt*(mcsel+reco+goodfit+livegate+rpitch+opa));
  ta->Project("+acc",binnames[ibin++],evtwt*(mcsel+reco+goodfit+livegate+rpitch+opa+rmom));

  double all = acc->GetBinContent(1);
  double norm = ngen;
  if(ngen < 0)
    norm = all;
  double prev = norm;
  for(ibin=1;ibin<=nbins;ibin++){
    racc->SetBinContent(ibin,acc->GetBinContent(ibin)/prev);
    prev = acc->GetBinContent(ibin);
  }
  cout << "Found " << norm << "Entries." << endl;
  racc->SetMaximum(1.1);
  acc->Scale(1.0/(float)norm);
  acc->SetMaximum(1.1);
  acc->SetStats(0);
  racc->SetStats(0);
  acc->GetXaxis()->SetLabelSize(0.06);
  racc->GetXaxis()->SetLabelSize(0.06);
  acc->SetMarkerSize(2.0);
  racc->SetMarkerSize(2.0);
  acc->GetYaxis()->SetTitleSize(0.05);
  racc->GetYaxis()->SetTitleSize(0.05);

  gStyle->SetPaintTextFormat("5.4f");
  TCanvas* acan = new TCanvas("acan","Acceptance",1200,800);
  acan->Clear();
  acan->Divide(1,2);
  acan->cd(1);
  TPad* tp = (TPad*)acan->cd(1);
  tp->SetBottomMargin(0.15);
  acc->Draw("histtext0");
  acan->cd(2);
  tp = (TPad*)acan->cd(2);
  tp->SetBottomMargin(0.15);
  racc->Draw("histtext0");
  if(strcmp(file,"")!=0)acan->SaveAs(file);
}

void hitres(TTree* ta) {
  TH1F* hresida = new TH1F("hresida","Hit Residual;mm",100,-2,2);
  TH1F* hresidna = new TH1F("hresidna","Hit Residual;mm",100,-2,2);
  TH1F* hresidall = new TH1F("hresidall","Hit Residual;mm",100,-2,2);
  hresida->SetFillColor(kGreen);
  hresidna->SetFillColor(kBlue);
  hresida->SetStats(0);
  hresidna->SetStats(0);
  THStack* hresidst = new THStack("hresidst","Hit Residual;mm");
  hresidst->Add(hresida);
  hresidst->Add(hresidna);
  ta->Project("hresida","detsh._resid","de.status>0&&detsh._active&&detsh._ambig!=0");
  ta->Project("hresidna","detsh._resid","de.status>0&&detsh._active&&detsh._ambig==0");
  ta->Project("hresidall","detsh._resid","de.status>0&&detsh._active");
  
  TH1F* hresa = new TH1F("hresa","Hit Drift Resolution;Reco R_{drift}-MC (mm)",100,-2,2);
  TH1F* hresna = new TH1F("hresna","Hit Drift Resolution;Reco R_{drift}-MC (mm)",100,-2,2);
  TH1F* hresall = new TH1F("hresall","Hit Drift Resolution;Reco R_{drift}-MC (mm)",100,-2,2);
  hresa->SetFillColor(kGreen);
  hresna->SetFillColor(kBlue);
  hresa->SetStats(0);
  hresna->SetStats(0);
  THStack* hresst = new THStack("hresst","Hit Drift Resolution;Reco R_{drift}-MC (mm)");
  hresst->Add(hresa);
  hresst->Add(hresna);
  ta->Project("hresa","detsh._rdrift-detshmc._dist","de.status>0&&detsh._active&&detsh._ambig!=0");
  ta->Project("hresna","detsh._rdrift-detshmc._dist","de.status>0&&detsh._active&&detsh._ambig==0");
  ta->Project("hresall","detsh._rdrift-detshmc._dist","de.status>0&&detsh._active");
  TCanvas* rescan = new TCanvas("rescan","rescan",800,800);
  rescan->Divide(1,2);
  rescan->cd(1);
  hresidst->Draw("h");
  hresidall->Fit("gaus","","sames");
  TLegend* rleg = new TLegend(0.15,0.6,0.4,0.85);
  rleg->AddEntry(hresida,"Resolved Ambiguity","f");
  rleg->AddEntry(hresidna,"Null Ambiguity","f");
  rleg->Draw();
  rescan->cd(2);
  hresst->Draw("h");
  hresall->Fit("gaus","","sames");
}

void wpull(TTree* ta) {
  TH1F* swp = new TH1F("swp","Final Fit Wire Position Pull",100,-25,25);
  TH1F* uwp = new TH1F("uwp","Final Fit Wire Position Pull",100,-25,25);
  TH1F* rwp = new TH1F("rwp","Final Fit Wire Position Pull",100,-25,25);
  swp->SetLineColor(kGreen);
  uwp->SetLineColor(kRed);
  rwp->SetLineColor(kBlue);
  ta->Project("swp","(detsh._wdist-detsh._hlen)/detsh._werr","de.status>0&&detsh._active&&detshmc._rel==0");
  ta->Project("uwp","(detsh._wdist-detsh._hlen)/detsh._werr","de.status>0&&detsh._active&&detshmc._rel<0");
  ta->Project("rwp","(detsh._wdist-detsh._hlen)/detsh._werr","de.status>0&&detsh._active&&detshmc._rel>0");
  TCanvas* wpcan = new TCanvas("wpcan","wpcan",800,800);
  wpcan->Divide(2,2);
  wpcan->cd(1);
  gPad->SetLogy();
  swp->Fit("gaus");
  wpcan->cd(2);
  uwp->Draw();
  wpcan->cd(3);
  gPad->SetLogy();
  rwp->Fit("gaus");
  TLegend* wleg = new TLegend(0.5,0.5,0.9,0.9);
  wleg->AddEntry(swp,"Primary Hit","L");
  wleg->AddEntry(uwp,"Unrelated Hit","L");
  wleg->AddEntry(rwp,"Related Hit","L");
  wpcan->cd(4);
  wleg->Draw();

}

void Ambig(TTree* ta,const char* file="") {
  gStyle->SetOptStat(1111);

  TCut ghit("detshmc._rel==0");
  TCut delta("detshmc._rel>0");
  TCut bkg("detshmc._rel<0");
  TCut gambig("detshmc._ambig==_ambig");
  TCut bambig("detshmc._ambig!=_ambig&&_ambig!=0");
  TCut nambig("detsh._ambig==0");
  TCut active("detsh._active>0");
// apply requested cuts

  TCut goodtrk("de.status>0&&de.nactive>20&&");


  TH1F* rdg = new TH1F("rdg","Hit fraction vs drift radius;true radius (mm);hit fraction",100,0.0,2.7);
  TH1F* rdn = new TH1F("rdn","Hit fraction vs drift radius;true radius (mm);hit fraction",100,0.0,2.7);
  TH1F* rdb = new TH1F("rdb","Hit fraction vs drift radius;true radius (mm);hit fraction",100,0.0,2.7);
  TH1F* rda = new TH1F("rda","Drift radius;true radius (mm)",100,0.0,2.7);
  TH1F* rdd = new TH1F("rdd","Drift radius;true radius (mm)",100,0.0,2.7);
  TH1F* rdf = new TH1F("rdf","Drift radius;true radius (mm)",100,0.0,2.7);
  TH1F* rdi = new TH1F("rdi","Drift radius;true radius (mm)",100,0.0,2.7);
  rdg->SetLineColor(kGreen);
  rdn->SetLineColor(kBlue);
  rdb->SetLineColor(kRed);
  rda->SetLineColor(kBlack);
  rdi->SetLineColor(kCyan);
  rdd->SetLineColor(kOrange);
  rdf->SetLineColor(kYellow);
  rdg->SetStats(0);
  rdn->SetStats(0);
//  rdb->SetStats(0);
  rdi->SetStats(0);
  rda->SetStats(0);
  rdd->SetStats(0);
  rdf->SetStats(0);
  rdg->Sumw2();
  rdn->Sumw2();
  rdb->Sumw2();
  rda->Sumw2();
  rdd->Sumw2();
  rdf->Sumw2();

  ta->Project("rdg","detshmc._dist",goodtrk+active+gambig+ghit);
  ta->Project("rdn","detshmc._dist",goodtrk+active+nambig+ghit);
  ta->Project("rdb","detshmc._dist",goodtrk+active+bambig+ghit);
  ta->Project("rda","detshmc._dist",goodtrk+active);
  ta->Project("rdi","detshmc._dist",goodtrk+ghit+(!active));
  ta->Project("rdd","detshmc._dist",goodtrk+active+delta);
  ta->Project("rdf","detshmc._dist",goodtrk+active+bkg);
  Double_t ntotal = rda->GetEntries();
  Double_t nright = rdg->GetEntries();
  Double_t nneutral = rdn->GetEntries();
  Double_t nwrong = rdb->GetEntries();
  Double_t ndelta = rdd->GetEntries();
  Double_t nbkg = rdf->GetEntries();
  Double_t ninact = rdi->GetEntries();
  TH1F* rdgr = new TH1F(*rdg);
  TH1F* rdnr = new TH1F(*rdn);
  TH1F* rdbr = new TH1F(*rdb);
  rdgr->Divide(rda);
  rdnr->Divide(rda);
  rdbr->Divide(rda);

  TH1F* momres0 = new TH1F("momres0","Momentum resolution at start of tracker;p_{reco}-p_{true}(MeV/c)",151,-4,4);
  TH1F* momres1 = new TH1F("momres1","Momentum resolution at start of tracker;p_{reco}-p_{true}(MeV/c)",151,-4,4);
  TH1F* momres2 = new TH1F("momres2","Momentum resolution at start of tracker;p_{reco}-p_{true}(MeV/c)",151,-4,4);
  momres0->SetLineColor(kBlack);
  momres1->SetMarkerColor(kCyan);
  momres1->SetMarkerStyle(4);
  momres2->SetMarkerColor(kOrange);
  momres2->SetMarkerStyle(5);
  ta->Project("momres0","de.mom-sqrt(demcent.momx^2+demcent.momy^2+demcent.momz^2)",goodtrk);
  ta->Project("momres1","de.mom-sqrt(demcent.momx^2+demcent.momy^2+demcent.momz^2)",goodtrk&&"de.status==1");
  ta->Project("momres2","de.mom-sqrt(demcent.momx^2+demcent.momy^2+demcent.momz^2)",goodtrk&&"de.status==2");

  TH1F* afg = new TH1F("afg","Average hit fraction vs momentum resolution;p_{reco}-p_{true}(MeV/c);hit fraction",41,-4,4);
  TH1F* afn = new TH1F("afn","Average hit fraction vs momentum resolution;p_{reco}-p_{true}(MeV/c):hit fraction",41,-4,4);
  TH1F* afb = new TH1F("afb","Average hit fraction vs momentum resolution;p_{reco}-p_{true}(MeV/c);hit fraction",41,-4,4);
  TH1F* afa = new TH1F("afa","Average hit fraction vs momentum resolution;p_{reco}-p_{true}(MeV/c);hit fraction",41,-4,4);
  afg->SetStats(0);
  afn->SetStats(0);
  afb->SetStats(0);
  afa->SetStats(0);
  afg->SetLineColor(kGreen);
  afn->SetLineColor(kBlue);
  afb->SetLineColor(kRed);
  afg->Sumw2();
  afn->Sumw2();
  afb->Sumw2();
  afa->Sumw2();
  ta->Project("afg","de.mom-sqrt(demcent.momx^2+demcent.momy^2+demcent.momz^2)",goodtrk+active+gambig);
  ta->Project("afn","de.mom-sqrt(demcent.momx^2+demcent.momy^2+demcent.momz^2)",goodtrk+active+nambig);
  ta->Project("afb","de.mom-sqrt(demcent.momx^2+demcent.momy^2+demcent.momz^2)",goodtrk+active+bambig);
  ta->Project("afa","de.mom-sqrt(demcent.momx^2+demcent.momy^2+demcent.momz^2)",goodtrk+active);
  afg->Divide(afa);
  afn->Divide(afa);
  afb->Divide(afa);
  afg->SetMinimum(0.0);
  afg->SetMaximum(1.1);

  TCanvas* ambigcan = new TCanvas("ambigcan","Hit Ambiguity",1200,800);
  ambigcan->Divide(2,2);


  ambigcan->cd(1);
  rda->Draw();
  rdi->Draw("same");
  rdd->Draw("same");
  rdf->Draw("same");
  TLegend* drleg = new TLegend(0.15,0.15,0.55,0.5);
  char dtitle[100];
  snprintf(dtitle,100,"%4.0f Active hits",ntotal);
  drleg->AddEntry(rda,dtitle,"l");
  snprintf(dtitle,100,"%4.4f Delta-ray hits",ndelta/ntotal);
  drleg->AddEntry(rdd,dtitle,"l");
  snprintf(dtitle,100,"%4.4f Inactive good hits",ninact/ntotal);
  drleg->AddEntry(rdi,dtitle,"l");
  snprintf(dtitle,100,"%4.4f Background hits",nbkg/ntotal);
  drleg->AddEntry(rdf,dtitle,"l");
  drleg->Draw();

  ambigcan->cd(2);

  rdgr->Draw();
  rdnr->Draw("same");
  rdbr->Draw("same");

  TLegend* leg = new TLegend(0.4,0.35,0.9,0.6);
  char title[80];
  snprintf(title,80,"Correct ambiguity %4.3f",nright/ntotal);
  leg->AddEntry(rdgr,title,"l");
  snprintf(title,80,"Null ambiguity %4.3f",nneutral/ntotal);
  leg->AddEntry(rdnr,title,"l");
  snprintf(title,80,"Incorrect ambiguity %4.3f",nwrong/ntotal);
  leg->AddEntry(rdbr,title,"l");
  leg->Draw();

  ambigcan->cd(3);
  double integral = momres0->GetEntries()*momres0->GetBinWidth(1);
  TF1* dscb = new TF1("dscb",fnc_dscb,-2.0,1.5,7);
  dscb->SetParName(0,"Norm");
  dscb->SetParName(1,"x0");
  dscb->SetParName(2,"sigma");
  dscb->SetParName(3,"ANeg");
  dscb->SetParName(4,"PNeg");
  dscb->SetParName(5,"APos");
  dscb->SetParName(6,"PPos");
  dscb->SetParameters(integral,0.0,0.15,1.0,4.5,1.2,10.0);
  dscb->SetNpx(1000);
  dscb->SetParLimits(3,0.0,50.0);
  dscb->SetParLimits(4,1.0,50.0);
  dscb->SetParLimits(5,0.0,50.0);
  dscb->SetParLimits(6,1.0,50.0);

  gPad->SetLogy();
  momres0->Fit("dscb","LIR");
  momres1->Draw("psame");
  momres2->Draw("psame");
  TLegend* mleg = new TLegend(0.13,0.6,0.43,0.85);
  snprintf(title,80,"%4.0f All Fits",momres0->GetEntries());
  mleg->AddEntry(momres0,title,"l");
  snprintf(title,80,"%4.0f Converged Fits",momres1->GetEntries());
  mleg->AddEntry(momres1,title,"p");
  snprintf(title,80,"%4.0f Unconverged Fits",momres2->GetEntries());
  mleg->AddEntry(momres2,title,"p");
  mleg->Draw();


  ambigcan->cd(4);
  afg->Draw();
  afn->Draw("same");
  afb->Draw("same");
  TLegend* fleg = new TLegend(0.16,0.35,0.625,0.6);
  fleg->AddEntry(afg,"Correct Ambiguity","l");
  fleg->AddEntry(afn,"No Ambiguity","l");
  fleg->AddEntry(afb,"Incorrect Ambiguity","l");
  fleg->Draw();

  ambigcan->cd(0);
  if(strcmp(file,"")!=0)ambigcan->SaveAs(file);
}

void Resid(TTree* ta) {

  TCut delta("detshmc._proc==17");
  TCut primary("detshmc._gen==2");
  TCut gambig("detshmc._ambig==detsh._ambig&&detsh._ambig!=0");
  TCut bambig("detshmc._ambig!=detsh._ambig&&detsh._ambig!=0");
  TCut nambig("detsh._ambig==0");
  TCut active("de.status==1 && detsh._active>0");
  TCut reco("de.status>0");
  TCut mcsel("sqrt(demcent.momx^2+demcent.momy^2+demcent.momz^2)>100.0");

  TH1F* rdg = new TH1F("rdg","True Drift radius;radius (mm);N hits",100,-0.05,2.55);
  TH1F* rdb = new TH1F("rdb","True Drift radius;radius (mm);N hits",100,-0.05,2.55);
  TH1F* rdn = new TH1F("rdn","True Drift radius;radius (mm);N hits",100,-0.05,2.55);
  rdg->SetLineColor(kBlue);
  rdb->SetLineColor(kRed);
  rdn->SetLineColor(kGreen);
  rdg->SetStats(0);
  rdb->SetStats(0);
  rdn->SetStats(0);

  TH1F* rpullg = new TH1F("rpullg","Correct Ambiguity Residual Pull;Pull;N hits",100,-6,6);
  TH1F* rpullb = new TH1F("rpullb","Incorrect Ambiguity Residual Pull;Pull;N hits",100,-6,6);
  TH1F* rpulln = new TH1F("rpulln","No Assigned Ambiguity Residual Pull;Pull;N hits",100,-6,6);
  TH1F* rpulld = new TH1F("rpulld","#delta-ray Residual Pull;Pull;N hits",100,-6,6);
  rpullg->SetLineColor(kBlue);
  rpullb->SetLineColor(kRed);
  rpulln->SetLineColor(kGreen);
  rpulld->SetLineColor(kCyan);

  ta->Project("rdg","detshmc._dist",reco+mcsel+active+gambig+primary);
  ta->Project("rdb","detshmc._dist",reco+mcsel+active+bambig+primary);
  ta->Project("rdn","detshmc._dist",reco+mcsel+active+nambig+primary);

  ta->Project("rpullg","detsh._resid/detsh._residerr",reco+mcsel+active+gambig+primary);
  ta->Project("rpullb","detsh._resid/detsh._residerr",reco+mcsel+active+bambig+primary);
  ta->Project("rpulln","detsh._resid/detsh._residerr",reco+mcsel+active+nambig+primary);
  ta->Project("rpulld","detsh._resid/detsh._residerr",reco+mcsel+active+delta);

  TCanvas* residcan = new TCanvas("residcan","Residuals",1200,800);
  residcan->Divide(2,1);
  residcan->cd(1);
  rdg->Draw();
  rdb->Draw("same");
  rdn->Draw("same");

  TLegend* leg = new TLegend(0.3,0.3,0.8,0.5);
  leg->AddEntry(rpullg,"Correct ambiguity","l");
  leg->AddEntry(rpullb,"Incorrect ambiguity","l");
  leg->AddEntry(rpulln,"No ambiguity assigned","l");
  leg->Draw();

  TPad* ppad = dynamic_cast<TPad*>(residcan->cd(2));
  ppad->Divide(1,4);
  ppad->cd(1);
  rpullg->Fit("gaus");
  ppad->cd(2);
  rpullb->Fit("gaus");
  ppad->cd(3);
  rpulln->Fit("gaus");
  ppad->cd(4);
  rpulld->Fit("gaus");

  residcan->cd(0);

}

void Con(TTree* ta) {
  TCut mcsel("sqrt(demcent.momx^2+demcent.momy^2+demcent.momz^2)>100.0");

  TH1F* con1 = new TH1F("con1","#chi^{2} fit consistency",500,0.0,1.0);
  TH1F* con2 = new TH1F("con2","#chi^{2} fit consistency",500,0.0,1.0);
  TH1F* lcon1 = new TH1F("lcon1","log_{10}(#chi^{2}) fit consistency",100,-10,0);
  TH1F* lcon2 = new TH1F("lcon2","log_{10}(#chi^{2}) fit consistency",100,-10,0);
  con1->SetLineColor(kBlue);
  con2->SetLineColor(kRed);
  lcon1->SetLineColor(kBlue);
  lcon2->SetLineColor(kRed);

  ta->Project("con1","de.fitcon",mcsel+"de.status==1");
  ta->Project("con2","de.fitcon",mcsel+"de.status==2");
  ta->Project("lcon1","log10(de.fitcon)",mcsel+"de.status==1");
  ta->Project("lcon2","log10(de.fitcon)",mcsel+"de.status==2");

  TCanvas* fcan = new TCanvas("fcan","fit consistency",500,800);
  fcan->Clear();
  fcan->Divide(1,2);
  fcan->cd(1);
  con1->Draw();
  con2->Draw("same");
  fcan->cd(2);
  lcon1->Draw();
  lcon2->Draw("same");

  TLegend* leg = new TLegend(0.1,0.6,0.4,0.8);
  leg->AddEntry(con1,"Fully Converged Fit","l");
  leg->AddEntry(con2,"Unconverged Fit","l");
  leg->Draw();

}

void TrkQualRes(TTree* ta,double tqcut) {
  TH1F* goodf = new TH1F("goodf","Momentum Resolution;Reco - True Momentum (MeV/c)",100,-4,4);
  TH1F* badf = new TH1F("badf","Momentum Resolution;Reco - True Momentum (MeV/c)",100,-4,4);
  goodf->SetLineColor(kGreen);
  goodf->SetLineWidth(2);
  goodf->SetFillColor(kGreen);
  goodf->SetFillStyle(3004);
  badf->SetLineColor(kMagenta);
  badf->SetFillColor(kMagenta);
  badf->SetLineWidth(2);
  badf->SetFillStyle(3005);
  goodf->SetStats(0);
  badf->SetStats(0);
  char tqcutgc[40];
  char tqcutbc[40];
  snprintf(tqcutgc,40,"detrkqual.trkqual>%f",tqcut);
  snprintf(tqcutbc,40,"detrkqual.trkqual<%f",tqcut);
  TCut tqcutg(tqcutgc);
  TCut tqcutb(tqcutbc);
  TCut reco("de.status>0");
  TCut mcsel("sqrt(demcent.momx^2+demcent.momy^2+demcent.momz^2)>100.0");
  ta->Project("goodf","de.mom-sqrt(demcent.momx^2+demcent.momy^2+demcent.momz^2)",reco+mcsel+tqcutg);
  ta->Project("badf","de.mom-sqrt(demcent.momx^2+demcent.momy^2+demcent.momz^2)",reco+mcsel+tqcutb);
  badf->Scale(goodf->GetEntries()/badf->GetEntries());
  TCanvas* tqcan = new TCanvas("tqcan","TrkQual",1000,800);
  goodf->Draw();
  badf->Draw("same");

  double gmax = goodf->GetMaximum();
  int gbin1 = goodf->FindFirstBinAbove(0.5*gmax);
  int gbin2 = goodf->FindLastBinAbove(0.5*gmax);
  double gfwhm = goodf->GetBinCenter(gbin2) - goodf->GetBinCenter(gbin1);

  double bmax = badf->GetMaximum();
  int bbin1 = badf->FindFirstBinAbove(0.5*bmax);
  int bbin2 = badf->FindLastBinAbove(0.5*bmax);
  double bfwhm = badf->GetBinCenter(bbin2) - badf->GetBinCenter(bbin1);

  char ltitle[40];

  TLegend* leg = new TLegend(0.55,0.6,0.9,0.9);
  snprintf(ltitle,40,"%s, FWHM =%3.0f KeV/c",tqcutgc,1000.0*gfwhm);
  leg->AddEntry(goodf,ltitle,"F");
  snprintf(ltitle,40,"%s, FWHM =%3.0f KeV/c",tqcutbc,1000.0*bfwhm);
  leg->AddEntry(badf,ltitle,"F");
  leg->Draw();
}

void StrawMat(TTree* ta) {
  TH1F* naddmat = new TH1F("naddmat","N Added Straws",30,-0.5,29.5);
  TH2F* matvshit = new TH2F("matvshit","N Straw vs N Hits",100,-0.5,99.5,100,-0.5,99.5);
  TH1F* addmatfrac = new TH1F("addmatfrac","Fraction of Added Straws",100,-0.01,0.5);
  addmatfrac->SetStats(0);
  matvshit->SetStats(0);

  TH1F* lofracres = new TH1F("lofracres","Momentum Resolution;Reco - True Mom. (MeV/c)",100,-2,2);
  TH1F* hifracres = new TH1F("hifracres","Momentum Resolution;Reco - True Mom. (MeV/c)",100,-2,2);
  lofracres->SetStats(0);
  hifracres->SetStats(0);
  lofracres->SetLineColor(kRed);
  hifracres->SetLineColor(kBlack);

  TH1F* hitdoca = new TH1F("hitdoca","DOCA to Wire;DOCA (mm)",100,-0.05,2.65);
  TH1F* adddoca = new TH1F("adddoca","DOCA to Wire;DOCA (mm)",100,-0.05,2.65);
  hitdoca->SetStats(0);
  hitdoca->SetLineColor(kRed);
  adddoca->SetStats(0);
  adddoca->SetLineColor(kBlue);

  TH1F* hitstraw = new TH1F("hitstraw","Straw Number;straw #",100,-0.5,99.5);
  TH1F* addstraw = new TH1F("addstraw","Straw number;straw #",100,-0.5,99.5);
  hitstraw->SetStats(0);
  hitstraw->SetLineColor(kRed);
  addstraw->SetStats(0);
  addstraw->SetLineColor(kBlue);

  ta->Project("addmatfrac","(de.nmatactive-de.nactive)/de.nmatactive","de.status>0");
  ta->Project("matvshit","de.nmatactive:de.nactive","de.status>0");
  ta->Project("naddmat","de.nmatactive-de.nactive","de.status>0");
  ta->Project("adddoca","detsm._doca","de.status>0&&detsm._active&&(!detsm._thita)");
  ta->Project("hitdoca","detsm._doca","de.status>0&&detsm._thita");
  ta->Project("addstraw","detsm._straw","de.status>0&&detsm._active&&(!detsm._thita)");
  ta->Project("hitstraw","detsm._straw","de.status>0&&detsm._thita");

  ta->Project("lofracres","de.mom-sqrt(demcent.momx^2+demcent.momy^2+demcent.momz^2)","de.status>0&&(de.nmatactive-de.nactive)/de.nmatactive<0.1");
  ta->Project("hifracres","de.mom-sqrt(demcent.momx^2+demcent.momy^2+demcent.momz^2)","de.status>0&&(de.nmatactive-de.nactive)/de.nmatactive>0.1");

  TLegend* leg = new TLegend(0.6,0.7,0.9,.9);
  leg->AddEntry(hitdoca,"Hit Straw","L");
  leg->AddEntry(adddoca,"Added Straw","L");

  TCanvas* mcan = new TCanvas("mcan","mcan",1200,800);
  mcan->Divide(3,2);
  mcan->cd(1);
  matvshit->Draw("colorz");
  mcan->cd(2);
  naddmat->Draw();
  mcan->cd(3);
  addmatfrac->Draw();
  mcan->cd(4);
  hitdoca->Draw();
  adddoca->Draw("same");
  leg->Draw();
  mcan->cd(5);
  hitstraw->Draw();
  addstraw->Draw("same");
  leg->Draw();
  mcan->cd(6);
  lofracres->Draw();
  hifracres->Draw("same");
  TLegend* mleg = new TLegend(0.6,0.7,0.9,0.9);
  mleg->AddEntry(lofracres,"N_{added}/N<0.1","L");
  mleg->AddEntry(hifracres,"N_{added}/N>0.1","L");
  mleg->Draw();
}

void TrkCaloHit(TTree* ta) {
  TCut goodtrkcalo("detrkqual.trkqual>0.6&&detch.active");
  TCut disk0("detch.disk==0");
  TCut disk1("detch.disk==1");
  TH1F* clen0 = new TH1F("clen0","TrkCaloHit POCA Crystal Depth;Depth(mm)",200,-50,250);
  TH1F* clen1 = new TH1F("clen1","TrkCaloHit POCA Crystal Depth;Depth(mm)",200,-50,250);
  TH1F* cdoca0 = new TH1F("cdoca0","TrkCaloHit DOCA;DOCA (mm)",200,-150,150);
  TH1F* cdoca1 = new TH1F("cdoca1","TrkCaloHit DOCA;DOCA (mm)",200,-150,150);
  TH1F* cdt0 = new TH1F("cdt0","TrkCaloHit #Delta t;TCH t_{0} - CaloCluster Time",100,-5,5);
  TH1F* cdt1 = new TH1F("cdt1","TrkCaloHit #Delta t;TCH t_{0} - CaloCluster Time",100,-5,5);
  TH1F* ep0 = new TH1F("ep0","TrkCaloHit E/P",100,0.0,1.1);
  TH1F* ep1 = new TH1F("ep1","TrkCaloHit E/P",100,0.0,1.1);
  TH1F* tdir0 = new TH1F("tdir0","Track Direction at POCA;#hat{t}#bullet#hat{#rho}",100,-1.0,1.0);
  TH1F* tdir1 = new TH1F("tdir1","Track Direction at POCA;#hat{t}#bullet#hat{#rho}",100,-1.0,1.0);
  TH1F* pr0 = new TH1F("pr0","POCA Radius;Transverse Radius (mm)",100,360,650);
  TH1F* pr1 = new TH1F("pr1","POCA Radius;Transverse Radius (mm)",100,360,650);

  clen0->SetLineColor(kRed);
  clen1->SetLineColor(kBlue);
  cdoca0->SetLineColor(kRed);
  cdoca1->SetLineColor(kBlue);
  cdt0->SetLineColor(kRed);
  cdt1->SetLineColor(kBlue);
  ep0->SetLineColor(kRed);
  ep1->SetLineColor(kBlue);
  tdir0->SetLineColor(kRed);
  tdir1->SetLineColor(kBlue);
  pr0->SetLineColor(kRed);
  pr1->SetLineColor(kBlue);

//  clen0->SetStats(0);
  clen1->SetStats(0);
//  cdoca0->SetStats(0);
  cdoca1->SetStats(0);
//  cdt0->SetStats(0);
  cdt1->SetStats(0);
//  ep0->SetStats(0);
  ep1->SetStats(0);
//  tdir0->SetStats(0);
  tdir1->SetStats(0);
//  pr0->SetStats(0);
  pr1->SetStats(0);
 
  ta->Project("clen0","detch.clen",goodtrkcalo&&disk0);
  ta->Project("clen1","detch.clen",goodtrkcalo&&disk1);
  ta->Project("cdoca0","detch.doca",goodtrkcalo&&disk0);
  ta->Project("cdoca1","detch.doca",goodtrkcalo&&disk1);
  ta->Project("cdt0","detch.t0-detch.ctime",goodtrkcalo&&disk0);
  ta->Project("cdt1","detch.t0-detch.ctime",goodtrkcalo&&disk1);
  ta->Project("ep0","detch.edep/sqrt(detch.momx^2+detch.momy^2+detch.momz^2)",goodtrkcalo&&disk0);
  ta->Project("ep1","detch.edep/sqrt(detch.momx^2+detch.momy^2+detch.momz^2)",goodtrkcalo&&disk1);
  ta->Project("tdir0","(detch.POCAx*detch.momx+detch.POCAy*detch.momy)/sqrt(detch.POCAx^2+detch.POCAy^2)/sqrt(detch.momx^2+detch.momy^2+detch.momz^2)",goodtrkcalo&&disk0);
  ta->Project("tdir1","(detch.POCAx*detch.momx+detch.POCAy*detch.momy)/sqrt(detch.POCAx^2+detch.POCAy^2)/sqrt(detch.momx^2+detch.momy^2+detch.momz^2)",goodtrkcalo&&disk1);
  ta->Project("pr0","sqrt(detch.POCAx^2+detch.POCAy^2)",goodtrkcalo&&disk0);
  ta->Project("pr1","sqrt(detch.POCAx^2+detch.POCAy^2)",goodtrkcalo&&disk1);

  
  TLegend* tchleg = new TLegend(0.6,0.7,0.9,0.9);
  tchleg->AddEntry(clen0,"Disk 0","L");
  tchleg->AddEntry(clen1,"Disk 1","L");

  TCanvas* tchcan = new TCanvas("tchcan","TrkCaloHit",800,600);
  tchcan->Divide(3,2);
  tchcan->cd(1);
  clen0->Draw();
  clen1->Draw("same");
  tchleg->Draw();
  tchcan->cd(2);
  cdoca0->Draw();
  cdoca1->Draw("same");
  tchcan->cd(3);
  cdt0->Draw();
  cdt1->Draw("same");
  tchcan->cd(4);
  ep0->Draw();
  ep1->Draw("same");
  tchcan->cd(5);
  tdir0->Draw();
  tdir1->Draw("same");
  tchcan->cd(6);
  pr0->Draw();
  pr1->Draw("same");
  

  TCut goodclen("detch.clen>0&&detch.clen<150.0");
  TCut badclen("detch.clen>150.0&&detch.clen<250.0");
//  TH1F* rad0g = new TH1F("rad0g","Cluster Radius, Disk 0",100,380,630);
//  TH1F* rad0b = new TH1F("rad0b","Cluster Radius, Disk 0",100,380,630);
//  ta->Project("rad0g","sqrt(detch.POCAx^2+detch.POCAy^2",goodtrkcalo&&disk0&&goodclen);
//  ta->Project("rad0b","sqrt(detch.POCAx^2+detch.POCAy^2",goodtrkcalo&&disk0&&badclen);
  TH2F* rad0 = new TH2F("rad0","POCA Radius vs Crystal Depth, Disk 0;Depth (mm);Radius (mm)",100,-50,250,100,360,650);
  TH2F* rad1 = new TH2F("rad1","POCA Radius vs Crystal Depth, Disk 1;Depth (mm);Radius (mm)",100,-50,250,100,360,650);
  TH2F* dot0 = new TH2F("dot0","Track Direction vs Crystal Depth, Disk 0;Depth (mm);#hat{t}#bullet#hat{#rho}",100,-50,250,100,-1,1);
  TH2F* dot1 = new TH2F("dot1","Track Direction vs Crystal Depth, Disk 1;Depth (mm);#hat{t}#bullet#hat{#rho}",100,-50,250,100,-1,1);
  TH2F* dotr0 = new TH2F("dotr0","Track Direction vs POCA Radius, Disk 0;Radius (mm);#hat{t}#bullet#hat{#rho}",100,360,650,100,-1,1);
  TH2F* dotr1 = new TH2F("dotr1","Track Direction vs POCA Radius, Disk 1;Radius (mm);#hat{t}#bullet#hat{#rho}",100,360,650,100,-1,1);
  TH2F* doca0 = new TH2F("doca0","DOCA vs Crystal Depth, Disk 0;Depth (mm);DOCA (mm)",100,-50,250,100,-100,100);
  TH2F* doca1 = new TH2F("doca1","DOCA vs Crystal Depth, Disk 1;Depth (mm);DOCA(mm)",100,-50,250,100,-100,100);
  TH2F* eopd0 = new TH2F("eopd0","E/P vs Crystal Depth, Disk 0;Depth (mm);E/P",100,-50,250,100,0.0,1.1);
  TH2F* eopd1 = new TH2F("eopd1","E/P vs Crystal Depth, Disk 1;Depth (mm);E/P",100,-50,250,100,0.0,1.1);
  TH2F* eopdir0 = new TH2F("eopdir0","E/P vs Track Direction, Disk 0;#hat{t}#bullet#hat{#rho};E/P",100,-1,1,100,0.0,1.1);
  TH2F* eopdir1 = new TH2F("eopdir1","E/P vs Track Direction, Disk 1;#hat{t}#bullet#hat{#rho};E/P",100,-1,1,100,0.0,1.1);
  TH2F* eopr0 = new TH2F("eopr0","E/P vs POCA Radius, Disk 0;Radius (mm);E/P",100,360,650,100,0.0,1.1);
  TH2F* eopr1 = new TH2F("eopr1","E/P vs POCA Radius, Disk 1;Radius (mm);E/P",100,360,650,100,0.0,1.1);
  rad0->SetStats(0);
  rad1->SetStats(0);
  dot0->SetStats(0);
  dot1->SetStats(0);
  eopd0->SetStats(0);
  eopd1->SetStats(0);
  eopdir0->SetStats(0);
  eopdir1->SetStats(0);
  eopr0->SetStats(0);
  eopr1->SetStats(0);
  ta->Project("rad0","sqrt(detch.POCAx^2+detch.POCAy^2):detch.clen",goodtrkcalo&&disk0);
  ta->Project("rad1","sqrt(detch.POCAx^2+detch.POCAy^2):detch.clen",goodtrkcalo&&disk1);
  ta->Project("dot0","(detch.POCAx*detch.momx+detch.POCAy*detch.momy)/sqrt(detch.POCAx^2+detch.POCAy^2)/sqrt(detch.momx^2+detch.momy^2+detch.momz^2):detch.clen",goodtrkcalo&&disk0);
  ta->Project("dot1","(detch.POCAx*detch.momx+detch.POCAy*detch.momy)/sqrt(detch.POCAx^2+detch.POCAy^2)/sqrt(detch.momx^2+detch.momy^2+detch.momz^2):detch.clen",goodtrkcalo&&disk1);
  ta->Project("dotr0","(detch.POCAx*detch.momx+detch.POCAy*detch.momy)/sqrt(detch.POCAx^2+detch.POCAy^2)/sqrt(detch.momx^2+detch.momy^2+detch.momz^2):sqrt(detch.POCAx^2+detch.POCAy^2)",goodtrkcalo&&disk0);
  ta->Project("dotr1","(detch.POCAx*detch.momx+detch.POCAy*detch.momy)/sqrt(detch.POCAx^2+detch.POCAy^2)/sqrt(detch.momx^2+detch.momy^2+detch.momz^2):sqrt(detch.POCAx^2+detch.POCAy^2)",goodtrkcalo&&disk1);
  ta->Project("eopd0","detch.edep/sqrt(detch.momx^2+detch.momy^2+detch.momz^2):detch.clen",goodtrkcalo&&disk0);
  ta->Project("eopd1","detch.edep/sqrt(detch.momx^2+detch.momy^2+detch.momz^2):detch.clen",goodtrkcalo&&disk1);
  ta->Project("eopdir0","detch.edep/sqrt(detch.momx^2+detch.momy^2+detch.momz^2):(detch.POCAx*detch.momx+detch.POCAy*detch.momy)/sqrt(detch.POCAx^2+detch.POCAy^2)/sqrt(detch.momx^2+detch.momy^2+detch.momz^2)",goodtrkcalo&&disk0);
  ta->Project("eopdir1","detch.edep/sqrt(detch.momx^2+detch.momy^2+detch.momz^2):(detch.POCAx*detch.momx+detch.POCAy*detch.momy)/sqrt(detch.POCAx^2+detch.POCAy^2)/sqrt(detch.momx^2+detch.momy^2+detch.momz^2)",goodtrkcalo&&disk1);
  ta->Project("eopr0","detch.edep/sqrt(detch.momx^2+detch.momy^2+detch.momz^2):sqrt(detch.POCAx^2+detch.POCAy^2)",goodtrkcalo&&disk0);
  ta->Project("eopr1","detch.edep/sqrt(detch.momx^2+detch.momy^2+detch.momz^2):sqrt(detch.POCAx^2+detch.POCAy^2)",goodtrkcalo&&disk1);
  ta->Project("doca0","detch.doca:detch.clen",goodtrkcalo&&disk0);
  ta->Project("doca1","detch.doca:detch.clen",goodtrkcalo&&disk1);
  TCanvas* tch0can = new TCanvas("tch0can","TrkCaloHit Disk 0", 1000, 800);
  tch0can->Divide(3,2);
  tch0can->cd(1);
  gPad->SetLogz();
  rad0->Draw("colorz");
  tch0can->cd(2);
  gPad->SetLogz();
  dot0->Draw("colorz");
  tch0can->cd(3);
  gPad->SetLogz();
  eopd0->Draw("colorz");
  tch0can->cd(4);
  gPad->SetLogz();
  doca0->Draw("colorz");
  tch0can->cd(5);
  gPad->SetLogz();
  eopdir0->Draw("colorz");
  tch0can->cd(6);
  gPad->SetLogz();
  eopr0->Draw("colorz");

  TCanvas* tch1can = new TCanvas("tch1can","TrkCaloHit Disk 1", 1000, 800);
  tch1can->Divide(3,2);
  tch1can->cd(1);
  gPad->SetLogz();
  rad1->Draw("colorz");
  tch1can->cd(2);
  gPad->SetLogz();
  dot1->Draw("colorz");
  tch1can->cd(3);
  gPad->SetLogz();
  eopd1->Draw("colorz");
  tch1can->cd(4);
  gPad->SetLogz();
  doca1->Draw("colorz");
  tch1can->cd(5);
  gPad->SetLogz();
  eopdir1->Draw("colorz");
  tch1can->cd(6);
  gPad->SetLogz();
  eopr1->Draw("colorz");

}

void TrkCaloHitMC(TTree* ta) {
  TCut goodtrkcalo("detrkqual.trkqual>0.6&&detch.active");
  TCut disk0("detch.disk==0");
  TCut disk1("detch.disk==1");

  TH1F* clen0m = new TH1F("clen0m","TrkCaloHit POCA Crystal Depth, Disk 0;Depth(mm)",200,-700,900);
  TH1F* clen0n = new TH1F("clen0n","TrkCaloHit POCA Crystal Depth, Disk 0;Depth(mm)",200,-700,900);
  TH1F* clen1m = new TH1F("clen1m","TrkCaloHit POCA Crystal Depth, Disk 1;Depth(mm)",200,-700,900);
  TH1F* clen1n = new TH1F("clen1n","TrkCaloHit POCA Crystal Depth, Disk 1;Depth(mm)",200,-700,900);
  TH1F* cdoca0m = new TH1F("cdoca0m","TrkCaloHit DOCA, Disk 0;DOCA (mm)",200,-250,700);
  TH1F* cdoca0n = new TH1F("cdoca0n","TrkCaloHit DOCA, Disk 0;DOCA (mm)",200,-250,700);
  TH1F* cdoca1m = new TH1F("cdoca1m","TrkCaloHit DOCA, Disk 1;DOCA (mm)",200,-250,700);
  TH1F* cdoca1n = new TH1F("cdoca1n","TrkCaloHit DOCA, Disk 1;DOCA (mm)",200,-250,700);
  clen0m->SetLineColor(kBlack);
  clen0n->SetLineColor(kCyan);
  clen1m->SetLineColor(kBlack);
  clen1n->SetLineColor(kCyan);
  cdoca0m->SetLineColor(kBlack);
  cdoca0n->SetLineColor(kCyan);
  cdoca1m->SetLineColor(kBlack);
  cdoca1n->SetLineColor(kCyan);
//
  clen0m->SetStats(0);
  clen0n->SetStats(0);
  clen1m->SetStats(0);
  clen1n->SetStats(0);
  cdoca0m->SetStats(0);
  cdoca0n->SetStats(0);
  cdoca1m->SetStats(0);
  cdoca1n->SetStats(0);

  TCut tcm("detchmc.prel>=0");
  TCut tcn("detchmc.prel<0");
  ta->Project("clen0m","detch.clen",goodtrkcalo&&disk0&&tcm);
  ta->Project("clen0n","detch.clen",goodtrkcalo&&disk0&&tcn);
  ta->Project("clen1m","detch.clen",goodtrkcalo&&disk1&&tcm);
  ta->Project("clen1n","detch.clen",goodtrkcalo&&disk1&&tcn);

  ta->Project("cdoca0m","detch.doca",goodtrkcalo&&disk0&&tcm);
  ta->Project("cdoca0n","detch.doca",goodtrkcalo&&disk0&&tcn);
  ta->Project("cdoca1m","detch.doca",goodtrkcalo&&disk1&&tcm);
  ta->Project("cdoca1n","detch.doca",goodtrkcalo&&disk1&&tcn);

  TLegend* tchmcleg = new TLegend(0.5,0.7,0.9,0.9);
  tchmcleg->AddEntry(clen0m,"Trk-Calo MC Match","L");
  tchmcleg->AddEntry(clen0n,"No MC Match","L");

  TCanvas* tchmc = new TCanvas("tchmc","TrkCaloHitMC",800,800);
  tchmc->Divide(2,2);
  tchmc->cd(1);
  gPad->SetLogy();
  clen0m->Draw();
  clen0n->Draw("same");
  tchmc->cd(2);
  gPad->SetLogy();
  clen1m->Draw();
  clen1n->Draw("same");
  tchmcleg->Draw();
  tchmc->cd(3);
  gPad->SetLogy();
  cdoca0m->Draw();
  cdoca0n->Draw("same");
  tchmc->cd(4);
  gPad->SetLogy();
  cdoca1m->Draw();
  cdoca1n->Draw("same");
}

void t0(TTree* ta) {
  TCut goodtrkcalo("detrkqual.trkqual>0.6&&detch.active");
  TCut goodtrknocalo("detrkqual.trkqual>0.6&&!detch.active");
  TCut disk0("detch.disk==0");
  TCut disk1("detch.disk==1");
  TH1F* t00 = new TH1F("t00","Track Fit t_{0} Resolution, TrkCaloHit;t_{0} reco - t_{0} MC (ns)",100,-5,5);
  TH1F* t01 = new TH1F("t01","Track Fit t_{0} Resolution, NoTrkCaloHit;t_{0} reco - t_{0} MC (ns)",100,-5,5);
//  t00->SetStats(0);
//  t01->SetStats(0);
//  t00->SetLineColor(kRed);
//  t01->SetLineColor(kBlue);
  ta->Project("t00","de.t0-fmod(demcmid.t0,1695)",goodtrkcalo);
  ta->Project("t01","de.t0-fmod(demcmid.t0,1695)",goodtrknocalo);
  TLegend* tchleg = new TLegend(0.6,0.7,0.9,0.9);
  tchleg->AddEntry(t00,"TrkCaloHit","L");
  tchleg->AddEntry(t01,"No TrkCaloHit","L");
  TCanvas* t0can = new TCanvas("t0can","t0can",600,400);
  t0can->Divide(2,1);
  t0can->cd(1);
  t00->Draw();
  t0can->cd(2);
  t01->Draw();
}

void Eff(TTree* ta, unsigned norm, double plo, double phi, int q=-1) {
  TH1F* allrec = new TH1F("allrec","Reco Fraction vs Generated Momentum",100,plo,phi);
  TH1F* tqrec = new TH1F("tqrec","Reco Fraction vs Generated Momentum",100,plo,phi);
  TH1F* t0rec = new TH1F("t0rec","Reco Fraction vs Generated Momentum",100,plo,phi);
  TH1F* tprtrig = new TH1F("tprtrig","Reco Fraction vs Generated Momentum",100,plo,phi);
  TH1F* cprtrig = new TH1F("cprtrig","Reco Fraction vs Generated Momentum",100,plo,phi);
  TH1F* trktrig = new TH1F("trktrig","Reco Fraction vs Generated Momentum",100,plo,phi);
  TH1F* alltrig = new TH1F("alltrig","Reco Fraction vs Generated Momentum",100,plo,phi);
  TH1F* cctrig = new TH1F("cctrig","Reco Fraction vs Generated Momentum",100,plo,phi);

  TCut t0cut("de.t0>700");
  TCut goodfit("detrkqual.trkqual>0.4");

  TCut cc("(trigbits&0x4)==0x4");
  // trigger depends on sign
  TCut goodtpr, goodcpr, goodtrk, goodtrig;
  if(q <0){
    // electrons
    goodtpr = TCut("(trigbits&0x200)==0x200");
    goodcpr = TCut("(trigbits&0x8)==0x8");
    goodtrk = TCut("(trigbits&0x208)>0");
    goodtrig = TCut("(trigbits&0x20C)>0");
  } else {
    // positrons
    goodtpr = TCut("(trigbits&0x400)==0x400");
    goodcpr = TCut("(trigbits&0x10)==0x10");
    goodtrk = TCut("(trigbits&0x410)>0");
    goodtrig = TCut("(trigbits&0x414)>0");
  }

  allrec->SetStats(0);
  tqrec->SetStats(0);
  tprtrig->SetStats(0);
  cprtrig->SetStats(0);
  trktrig->SetStats(0);
  alltrig->SetStats(0);
  cctrig->SetStats(0);

  allrec->SetLineColor(kBlue);
  tqrec->SetLineColor(kRed);
  tprtrig->SetLineColor(kGreen);
  cprtrig->SetLineColor(kOrange);
  trktrig->SetLineColor(kCyan);
  alltrig->SetLineColor(kBlack);
  cctrig->SetLineColor(kYellow);
  ta->Project("allrec","sqrt(demcgen.momx^2+demcgen.momy^2+demcgen.momz^2)");
  ta->Project("tqrec","sqrt(demcgen.momx^2+demcgen.momy^2+demcgen.momz^2)",goodfit);
  ta->Project("tprtrig","sqrt(demcgen.momx^2+demcgen.momy^2+demcgen.momz^2)",goodfit&&goodtpr);
  ta->Project("cprtrig","sqrt(demcgen.momx^2+demcgen.momy^2+demcgen.momz^2)",goodfit&&goodcpr);
  ta->Project("trktrig","sqrt(demcgen.momx^2+demcgen.momy^2+demcgen.momz^2)",goodfit&&goodtrk);
  ta->Project("alltrig","sqrt(demcgen.momx^2+demcgen.momy^2+demcgen.momz^2)",goodfit&&goodtrig);
  ta->Project("cctrig","sqrt(demcgen.momx^2+demcgen.momy^2+demcgen.momz^2)",goodfit&&cc);

  // scale by the absolute normalization
  double scalefac =100.0/norm;
  allrec->Scale(scalefac);
  tqrec->Scale(scalefac);
  tprtrig->Scale(scalefac);
  cprtrig->Scale(scalefac);
  trktrig->Scale(scalefac);
  alltrig->Scale(scalefac);
  cctrig->Scale(scalefac);
  TCanvas* effcan = new TCanvas("effcan","Efficiency",600,600);
  allrec->Draw("h");
  tqrec->Draw("hsame");
  alltrig->Draw("hsame");
  trktrig->Draw("hsame");
  tprtrig->Draw("hsame");
  cprtrig->Draw("hsame");
  alltrig->Draw("hsame");
  cctrig->Draw("hsame");
  TLegend* leg = new TLegend(0.1,0.6,0.5,0.9);
  leg->AddEntry(allrec,"All Reco","l");
  leg->AddEntry(tqrec,"TrkQual>0.4","l");
  leg->AddEntry(alltrig,"All Trigger","l");
  leg->AddEntry(trktrig,"Track Trigger","l");
  leg->AddEntry(tprtrig,"TrackPatRec Trigger","l");
  leg->AddEntry(cprtrig,"CalPatRec Trigger","l");
  leg->AddEntry(cctrig,"CaloCluster Trigger","l");
  leg->Draw();
}

void PlotIPA(TTree* ta) {
  gStyle->SetOptFit(111111);
  gStyle->SetOptStat("oumr");
  TH1F* trkqual = new TH1F("trkqual","TrkQual",103,-0.01,1.01);
  TH1F* mom = new TH1F("mom","Reco Momentum;P_{reco} (MeV/c)",100,40,56);
  TH1F* momres = new TH1F("momres","Momentum Resolution;P_{reco} - P_{MC} (MeV/c)",100,-5.0,5.0);
  TH1F* nactive = new TH1F("nactive","N Active Straw Hits",121,-0.5,120.5);
  trkqual->SetStats(0);
  ta->Project("trkqual","de.trkqual");
  ta->Project("mom","de.mom","de.trkqual>0.4");
  ta->Project("momres","de.mom-sqrt(demcent.momx^2+demcent.momy^2+demcent.momz^2)","de.trkqual>0.4");
  ta->Project("nactive","de.nactive","de.trkqual>0.4");
  TCanvas* ipacan = new TCanvas("ipacan","ipacan",800,800);
  ipacan->Divide(2,2);
  ipacan->cd(1);
  trkqual->Draw();
  ipacan->cd(2);
  nactive->Draw();
  ipacan->cd(3);
  mom->Draw();
  ipacan->cd(4);
  gPad->SetLogy();
  double integral = momres->GetEntries()*momres->GetBinWidth(1);
  cout << "Integral = " << integral << " mean = " << momres->GetMean() << " rms = " << momres->GetRMS() << endl;
  TF1* dscb = new TF1("dscb",fnc_dscb,-10.0,5,7);
  dscb->SetParName(0,"Norm");
  dscb->SetParName(1,"x0");
  dscb->SetParName(2,"sigma");
  dscb->SetParName(3,"ANeg");
  dscb->SetParName(4,"PNeg");
  dscb->SetParName(5,"APos");
  dscb->SetParName(6,"PPos");
  dscb->SetParameters(2*integral,0.0,0.15,1.0,4.5,1.2,10.0);
  momres->Fit("dscb","RQ");
}
