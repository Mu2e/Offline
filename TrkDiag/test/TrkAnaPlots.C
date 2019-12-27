#include "TH1F.h"
#include "TF1.h"
#include "TFile.h"
#include "TDirectory.h"
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

class TrkAnaPlots {
  public:
  enum TrackerRegion {entrance=0,middle, exit};
  TrkAnaPlots(TFile* file,float momwin=1.5) {
    _tn = (TTree*)( ((TDirectory*)file->Get("TrkAnaNeg"))->Get("trkana"));
    _tp = (TTree*)( ((TDirectory*)file->Get("TrkAnaPos"))->Get("trkana"));
    BuildCuts(momwin);
    _rname.push_back("Entrance");
    _rname.push_back("Middle");
    _rname.push_back("Exit");
  }
  
  void BuildCuts(float momwin);
  void PID();
  void FitMomResp(TH1F* momresp);
  void FitMomRes(TH1F* momres);
  void MomResp(TrackerRegion region=entrance);
  void MomRes(TrackerRegion region=entrance);
  void SelPlots(int charge=-1);
  void Acc(int ngen,int charge=-1);
  void hitres();
  void wpull();
  void Ambig();
  void Resid();
  void Con();
  void TrkQual(const char* extra="");
  void TrkQualRes(float tqcut);
  void StrawMat();
  void TrkCaloHit(float tqcut,int pdg=11);
  void TrkCaloHitMC();
  void dEdx();
  void t0();
  void Eff(unsigned norm, double plo, double phi, int q=-1);
  void PlotIPA();
  void Upstream();
  void PBI(unsigned ngen, int charge=-1);
  void Trigger();
  void Alg();
// cuts
  TCut _reco, _goodfit, _rpitch, _livegate, _opa, _upstream, _physics, _final, _pbi;
  TCut _eminus,_eplus,_ele, _muminus, _muplus, _mu;
  TCut _CRV, _eminustrig, _eplustrig, _eminusrmom, _eplusrmom, _eminuspid, _epluspid, _eminustq, _eplustq, _downstream;
  TCut _TPR, _CPR;
 // Trees 
  TTree* _tn;
  TTree* _tp; 
  // canvases
  TCanvas *_pidcan, *_pidqcan, *_pidmomcan, *_radcan, *_rscan, *_rcan, *_acan, *_ecan, *_rescan, *_wpcan, *_ambigcan, *_residcan, *_fcan, *_tqcan, *_tqrcan, *_mcan, *_tchcan, *_tch0can, *_tch1can, *_dtchtcan,*_t0can, *_effcan, *_ipacan, *_ucan, *_uecan, *_tchmccan, *_spcan;
  vector<string> _rname;
};

void TrkAnaPlots::BuildCuts(float momwin){
  char ctext[200];

  _reco = TCut("de.status>0");
  _ele = TCut("abs(demc.pdg)==11");
  _eminus = TCut("demc.pdg==11");
  _eplus = TCut("demc.pdg==-11");
  _mu = TCut("abs(demc.pdg)==13");
  _muminus = TCut("demc.pdg==13");
  _muplus = TCut("demc.pdg==-13");
  double tdlow(0.57735027);
  double tdhigh(1.0);
  double t0min(700.0);
  double t0max(1650.0);
  double eminusmom(105.0);
  double eplusmom(92.3);
  snprintf(ctext,200,"deent.td>%5.5f&&deent.td<%5.5f",tdlow,tdhigh);
  _rpitch = TCut(ctext);
  snprintf(ctext,200,"de.t0>%f&&de.t0<%f",t0min,t0max);
  _livegate = TCut(ctext);
  _opa = TCut("abs(deent.d0<105) && abs(deent.d0+2/deent.om)>450 && abs(deent.d0+2/deent.om)<680");
  _upstream = TCut("ue.status<0");
  _pbi = TCut("evtwt.PBIWeight");
  double crvlow(-50.0);
  double crvhigh(150.0);
  snprintf(ctext,200,"bestcrv<0||(de.t0-crvinfo._timeWindowStart[bestcrv]<%5.1f||de.t0-crvinfo._timeWindowStart[bestcrv]>%5.1f)",crvlow,crvhigh);
  _CRV = TCut(ctext);
  _eminustq = TCut("dequal.TrkQualDeM>0.8");
  _eplustq = TCut("dequal.TrkQualDeP>0.8");
  _eminuspid = TCut("dequal.TrkPIDDeM>0.95");
  _epluspid = TCut("dequal.TrkPIDDeP>0.95");
  _eminustrig = TCut("(trigbits&0x4010)>0"); // negative, TrkPatRec or CalPatRec
  _eplustrig = TCut("(trigbits&0x8020)>0"); // positive, TrkPatRec or CalPatRec
  snprintf(ctext,200,"abs(deent.mom-%f)<%f",eminusmom,momwin);
  _eminusrmom = TCut(ctext);
  snprintf(ctext,200,"abs(deent.mom-%f)<%f",eplusmom,momwin);
  _eplusrmom = TCut(ctext);
  _downstream = TCut("demcent.momz>0");
  _physics = _rpitch+_opa+_livegate;
  _TPR = TCut("de.alg==0");
  _CPR = TCut("de.alg==1");
  _final = _eminustrig+_reco+_eminustq+_livegate+_rpitch+_opa+_upstream+_CRV+_eminuspid+_eminusrmom;
}

void TrkAnaPlots::PID() {
  TH2F* fr1vr0nc = new TH2F("fr1vr0nc","Disk 1 front projected radius vs Disk 0 front projected radius, no TrkCaloHit;Front R_{0} (mm);Front R_{1} (mm)",100,0.0,700.0,100,0.0,700);
  TH2F* fr1vr0c = new TH2F("fr1vr0c","Disk 1 front projected radius vs Disk 0 front projected radius, TrkCaloHit;Front R_{0} (mm);Front R_{1} (mm)",100,0.0,700.0,100,0.0,700);
  TH2F* br1vr0nc = new TH2F("br1vr0nc","Disk 1 back projected radius vs Disk 0 back projected radius, no TrkCaloHit;Back R_{0} (mm);Back R_{1} (mm)",100,0.0,700.0,100,0.0,700);
  TH2F* br1vr0c = new TH2F("br1vr0c","Disk 1 back projected radius vs Disk 0 back projected radius, TrkCaloHit;Back R_{0} (mm);Back R_{1} (mm)",100,0.0,700.0,100,0.0,700);
  TH1F* dempid = new TH1F("dempid","TCHPidQual",120,-0.1,1.1);
  TH1F* deppid = new TH1F("deppid","TCHPidQual",120,-0.1,1.1);
  TH1F* dmumpid = new TH1F("dmumpid","TCHPidQual",120,-0.1,1.1);
  TH1F* dmuppid = new TH1F("dmuppid","TCHPidQual",120,-0.1,1.1);
  TH1F* demom = new TH1F("demom","Momentum;Mom (MeV/c)",120,60,180);
  TH1F* dmmom = new TH1F("dmmom","Momentum;Mom (MeV/c)",120,60,180);
  TH1F* demompid = new TH1F("demompid","Momentum after PIDQual cut;Mom (MeV/c)",120,60,180);
  TH1F* dmmompid = new TH1F("dmmompid","Momentum after PIDQual cut;Mom (MeV/c)",120,60,180);
  TH1F* demompide = new TH1F("demompide","PIDQual efficiency vs Momentum;Mom (MeV/c)",120,60,180);
  TH1F* dmmompide = new TH1F("dmmompide","PIDQual efficiency vs Momentum;Mom (MeV/c)",120,60,180);

  TH2F* deevsp = new TH2F("deevsp","E_{calo} vs P_{track};P (MeV/c);E (MeV)",50,50,200,50,0,220);
  TH2F* dmevsp = new TH2F("dmevsp","E_{calo} vs P_{track};P (MeV/c);E (MeV)",50,50,200,50,0,220);
  TH2F* dpevsp = new TH2F("dpevsp","E_{calo} vs P_{track};P (MeV/c);E (MeV)",50,50,200,50,0,220);
  TProfile* deevspp = new TProfile("deevspp","E_{calo} vs P_{track};P (MeV/c);E (MeV)",50,50,200,0,220);
  TProfile* dmevspp = new TProfile("dmevspp","E_{calo} vs P_{track};P (MeV/c);E (MeV)",50,50,200,0,220);
  TProfile* dpevspp = new TProfile("dpevspp","E_{calo} vs P_{track};P (MeV/c);E (MeV)",50,50,200,0,220);


  fr1vr0nc->SetStats(0);
  fr1vr0c->SetStats(0);
  br1vr0nc->SetStats(0);
  br1vr0c->SetStats(0);
  dempid->SetStats(0);
  deppid->SetStats(0);
  dmumpid->SetStats(0);
  dmuppid->SetStats(0);
  //
  dempid->SetLineColor(kBlue);
  deppid->SetLineColor(kGreen);
  dmumpid->SetLineColor(kMagenta);
  dmuppid->SetLineColor(kRed);

  demom->Sumw2();
  dmmom->Sumw2();
  demompid->Sumw2();
  dmmompid->Sumw2();
  demompide->Sumw2();
  dmmompide->Sumw2();
  demom->SetLineColor(kRed);
  dmmom->SetLineColor(kBlue);
  demompid->SetLineColor(kRed);
  dmmompid->SetLineColor(kBlue);
  demompide->SetLineColor(kRed);
  dmmompide->SetLineColor(kBlue);
  deevsp->SetStats(0);
  dmevsp->SetStats(0);
  deevspp->SetStats(0);
  dmevspp->SetStats(0);
  deevsp->SetLineColor(kRed);
  dmevsp->SetLineColor(kBlue);
  deevspp->SetLineColor(kRed);
  dmevspp->SetLineColor(kBlue);
  
  demompide->SetStats(0);
  dmmompide->SetStats(0);
  
  _tn->Project("dempid","detrkpid.mvaout",_downstream&&_eminus);
  _tn->Project("dmumpid","detrkpid.mvaout",_downstream&&_muminus);
  _tp->Project("deppid","detrkpid.mvaout",_downstream&&_eplus);
  _tp->Project("dmuppid","detrkpid.mvaout",_downstream&&_muplus);
  _tn->Project("deevsp","detch.edep:dexit.mom",_downstream&&_eminus);
  _tp->Project("+deevsp","detch.edep:dexit.mom",_downstream&&_eplus);
  _tn->Project("deevspp","detch.edep:dexit.mom",_downstream&&_eminus);
  _tp->Project("+deevspp","detch.edep:dexit.mom",_downstream&&_eplus);
  _tn->Project("dmevsp","detch.edep:dexit.mom",_downstream&&_muminus);
  _tp->Project("+dmevsp","detch.edep:dexit.mom",_downstream&&_muplus);
  _tn->Project("dmevspp","detch.edep:dexit.mom",_downstream&&_muminus);
  _tp->Project("+dmevspp","detch.edep:dexit.mom",_downstream&&_muplus);
  _tn->Project("dpevsp","detch.edep:dexit.mom",_upstream&&_eplus);
  _tp->Project("+dpevsp","detch.edep:dexit.mom",_upstream&&_eminus);
  _tn->Project("dpevspp","detch.edep:dexit.mom",_upstream&&_eplus);
  _tp->Project("+dpevspp","detch.edep:dexit.mom",_upstream&&_eminus);

  _tn->Project("fr1vr0c","detrkpid.disk1frad:detrkpid.disk0frad",_downstream&&_eminus&&"detch.active");
  _tn->Project("fr1vr0nc","detrkpid.disk1frad:detrkpid.disk0frad",_downstream&&_eminus&&"!detch.active");
  _tp->Project("+fr1vr0c","detrkpid.disk1frad:detrkpid.disk0frad",_downstream&&_eplus&&"detch.active");
  _tp->Project("+fr1vr0nc","detrkpid.disk1frad:detrkpid.disk0frad",_downstream&&_eplus&&"!detch.active");

  _tn->Project("br1vr0c","detrkpid.disk1brad:detrkpid.disk0brad",_downstream&&_eminus&&"detch.active");
  _tn->Project("br1vr0nc","detrkpid.disk1brad:detrkpid.disk0brad",_downstream&&_eminus&&"!detch.active");
  _tp->Project("+br1vr0c","detrkpid.disk1brad:detrkpid.disk0brad",_downstream&&_eplus&&"detch.active");
  _tp->Project("+br1vr0nc","detrkpid.disk1brad:detrkpid.disk0brad",_downstream&&_eplus&&"!detch.active");

  _pidcan = new TCanvas("pidcan","pidcan",800,1200);
  _pidcan->Divide(1,2);
  _pidcan->cd(1);
  gPad->SetLogy();
  if(dempid->GetEntries() > deppid->GetEntries()){
    dempid->Draw();
    deppid->Draw("same");
  } else {
    deppid->Draw();
    dempid->Draw("same");
  }
  dmumpid->Draw("same");
  dmuppid->Draw("same");
  
  TLegend* leg = new TLegend(0.1,0.6,0.5,0.9);
  leg->AddEntry(dempid,"True e^{-}","L");
  leg->AddEntry(dmumpid,"True #mu^{-}","L");
  leg->AddEntry(deppid,"True e^{+}","L");
  leg->AddEntry(dmuppid,"True #mu^{+}","L");
  leg->Draw();

  _radcan = new TCanvas("radcan","radcan",800,800);
  _radcan->Divide(2,2);
  _radcan->cd(1);
  fr1vr0c->Draw("colorz");
  _radcan->cd(2);
  fr1vr0nc->Draw("colorz");
  _radcan->cd(3);
  br1vr0c->Draw("colorz");
  _radcan->cd(4);
  br1vr0nc->Draw("colorz");

  TH2F* qvqmu = new TH2F("qvqmu","TrkPID vs TrkQual, true #mu^{#pm};TrkQual mvaout;TrkPid mvaout",50,-0.05,1.05,50,-0.05,1.05);
  TH2F* qvqe = new TH2F("qvqe","TrkPID vs TrkQual, true e^{#pm};TrkQual mvaout;TrkPid mvaout",50,-0.05,1.05,50,-0.05,1.05);
  qvqmu->SetStats(0);
  qvqe->SetStats(0);

  _tn->Project("qvqmu","dequal.TrkQualDeM:dequal.TrkPIDDeM",_downstream&&_muminus);
  _tp->Project("+qvqmu","dequal.TrkQualDeP:dequal.TrkPIDDeP",_downstream&&_muplus);
  _tn->Project("qvqe","dequal.TrkQualDeM:dequal.TrkPIDDeM",_downstream&&_eminus);
  _tp->Project("+qvqe","dequal.TrkQualDeP:dequal.TrkPIDDeP",_downstream&&_eplus);

  _pidqcan = new TCanvas("pidqcan","pidqcan",1000,500);
  _pidqcan->Divide(2,1);
  _pidqcan->cd(1);
  gPad->SetLogz();
  qvqe->Draw("colorz");
  _pidqcan->cd(2);
  gPad->SetLogz();
  qvqmu->Draw("colorz");

  
  _tn->Project("dmmom","deent.mom",_downstream&&_muminus&&_eminustq);
  _tp->Project("+dmmom","deent.mom",_downstream&&_muplus&&_eplustq);
  _tn->Project("demom","deent.mom",_downstream&&_eminus&&_eminustq);
  _tp->Project("+demom","deent.mom",_downstream&&_eplus&&_eplustq);
  _tn->Project("dmmompid","deent.mom",_downstream&&_muminus&&_eminustq&&_eminuspid);
  _tp->Project("+dmmompid","deent.mom",_downstream&&_muplus&&_eplustq&&_epluspid);
  _tn->Project("demompid","deent.mom",_downstream&&_eminus&&_eminustq&&_eminuspid);
  _tp->Project("+demompid","deent.mom",_downstream&&_eplus&&_eplustq&&_epluspid);

  demompide->Divide(demompid,demom);
  dmmompide->Divide(dmmompid,dmmom);

  _pidmomcan = new TCanvas("mompid","mompid",800,800);
  _pidmomcan->Divide(2,2);
  _pidmomcan->cd(1);
  demom->Draw();
  dmmom->Draw("same");
  TLegend* mpleg = new TLegend(0.1,0.7,0.4,0.9);
  mpleg->AddEntry("demom","True electron","L");
  mpleg->AddEntry("dmmom","True muon","L");
  mpleg->Draw();
  _pidmomcan->cd(2);
  demompid->Draw();
  dmmompid->Draw("same");
  _pidmomcan->cd(3);
  demompide->SetMinimum(-0.05);
  demompide->SetMaximum(1.05);
  demompide->Draw();
  dmmompide->Draw("same");
  _pidmomcan->cd(4);
  deevsp->Draw("box");
  deevspp->Fit("pol1","","same");
  dmevsp->Draw("boxsame");
  dmevspp->Fit("pol1","","same");
  mpleg->Draw();
}

void TrkAnaPlots::FitMomResp(TH1F* momresp) {
  double integral = momresp->GetEntries()*momresp->GetBinWidth(1);
  cout << "Integral = " << integral << " mean = " << momresp->GetMean() << " rms = " << momresp->GetRMS() << endl;
  TF1* dscb = new TF1("dscb",fnc_dscb,-10.0,5,7);
  dscb->SetParName(0,"Norm");
  dscb->SetParName(1,"x0");
  dscb->SetParName(2,"sigma");
  dscb->SetParName(3,"ANeg");
  dscb->SetParName(4,"PNeg");
  dscb->SetParName(5,"APos");
  dscb->SetParName(6,"PPos");

  dscb->SetParameters(integral,momresp->GetMean()+0.5,0.3,0.7,3.0,3.0,6.0);
  dscb->SetNpx(1000);
  dscb->SetParLimits(2,0.1,50.0);
  dscb->SetParLimits(3,0.0,50.0);
  dscb->SetParLimits(4,1.0,50.0);
  dscb->SetParLimits(5,0.0,50.0);
  dscb->SetParLimits(6,1.0,50.0);

  momresp->SetMinimum(0.5);
  momresp->Fit("dscb","LRQ");
  momresp->Fit("dscb","LRM");
}

void TrkAnaPlots::MomResp(TrackerRegion region) {
  gStyle->SetOptFit(111111);
  gStyle->SetOptStat("oumr");
  char title[80];
  snprintf(title,80,"Momentum Response at Tracker %s;P_{reconstructed} - P_{generated} (MeV/c)",_rname[region].c_str());
  TH1F* momresp = new TH1F("momresp",title,151,-4.0,4.0);
  TH1F* momrespr = new TH1F("momrespr",title,151,-4.0,4.0);
  momresp->Sumw2();
  momrespr->Sumw2();
  momresp->SetLineColor(kRed);
  momrespr->SetLineColor(kBlue);
  _tn->Project("momresp","deent.mom-sqrt(demcgen.momx^2+demcgen.momy^2+demcgen.momz^2)",_final*_pbi);
  _tn->Project("momrespr","deent.mom-sqrt(demcgen.momx^2+demcgen.momy^2+demcgen.momz^2)",(_physics+_eminusrmom)*_pbi);
  _rscan = new TCanvas("rscan","Momentum Response",800,800);
  gPad->SetLogy();
  FitMomResp(momrespr);
  FitMomResp(momresp);
  momrespr->GetFunction("dscb")->SetLineColor(kBlue);
  momrespr->Draw();
  momrespr->SetStats(0);
  momresp->Draw("same");
  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9);
  leg->AddEntry(momrespr,"All Tracks","l");
  leg->AddEntry(momresp,"Selected Tracks","l");
  leg->Draw();
}

void TrkAnaPlots::FitMomRes(TH1F* momres) {
  TF1* dscb = new TF1("dscb",fnc_dscb,-2.0,2.5,7);
  dscb->SetParName(0,"Norm");
  dscb->SetParName(1,"x0");
  dscb->SetParName(2,"sigma");
  dscb->SetParName(3,"ANeg");
  dscb->SetParName(4,"PNeg");
  dscb->SetParName(5,"APos");
  dscb->SetParName(6,"PPos");

  double integral = momres->GetEntries()*momres->GetBinWidth(1);
  cout << "Integral = " << integral << " mean = " << momres->GetMean() << " rms = " << momres->GetRMS() << endl;
  dscb->SetParameters(2*integral,0.0,0.15,1.0,4.5,1.2,10.0);

  momres->Fit("dscb","RQ");
  momres->Fit("dscb","RMQ");
  TFitResultPtr fitres = momres->Fit("dscb","SRE");
  cout << "Core Sigma = " << dscb->GetParameter(2) << " +- " << dscb->GetParError(2) << endl;
  cout << "High Side Power = " << dscb->GetParameter(6) << " +- " << dscb->GetParError(6) << endl;

  momres->Fit("dscb","RQ");
  momres->Fit("dscb","RMQ");

  // count outliers
  int outbin = momres->FindBin(1.1);
  double outint = momres->Integral(outbin,momres->GetNbinsX());
  double totint = momres->Integral();
  double outrat = outint/totint;
  cout <<"Outlier integral = " << outint  << " total integral = " << totint << " Outlier fraction = " << outrat << endl;

  TLine* zero = new TLine(0.0,0.0,0.0,momres->GetBinContent(momres->GetMaximumBin()));
  zero->SetLineStyle(2);
  zero->Draw();

//  TPaveText* rtext = new TPaveText(0.1,0.5,0.4,0.9,"NDC");
//  rtext->AddText("Reco Cuts");
//  char line[80];
//  snprintf(line,80,"%s",cut.GetTitle());
//  rtext->AddText(line);
//  sprintf(line,"%5.0f Tracks",momres->GetEntries());
//  rtext->AddText(line);
//  rtext->Draw();
}

void TrkAnaPlots::MomRes(TrackerRegion region) {
// cuts
  char title[80];
  snprintf(title,80,"momentum resolution at tracker %s;MeV/c",_rname[region].c_str());
  TH1F* momres = new TH1F("momres",title,251,-4,4);
  momres->Sumw2();
  _tn->Project("momres","deent.mom-sqrt(demcent.momx^2+demcent.momy^2+demcent.momz^2)",_final*_pbi);
//  momres->SetMinimum(0.5);
  _rcan = new TCanvas("rcan","Momentum Resolution",800,800);
  gStyle->SetOptFit(111111);
  gStyle->SetOptStat("oumr");
  gPad->SetLogy();
  FitMomRes(momres);

}

void TrkAnaPlots::SelPlots(int charge) {
  vector<string> names={"trig", "trkqual", "t0", "pitch","d0", "rmax", "PID", "CRVDT", "Momentum"};
  vector<string> titles={"Trigger Bits;log2(trigger)", "Track Fit Quality;trkqual", "t0", "Reco pitch;tan(#lambda)", "d0;d0 (mm)", "Rmax;rmax (mm)", "PID;TrkPIDDeM", "CRV #Deltat;t0-#t_{CRV}", "Track Momentum;Momentum (MeV/c)"};
  vector<string> vars={"log2(trigbits)", "dequal.TrkQualDeM", "de.t0", "deent.td", "deent.d0", "abs(de.ent.d0+2.0/deent.om)", "dequal.TrkPIDDeM", "de.t0-crvinfo._timeWindowStart[bestcrv]", "deent.mom"};
  vector<double> low={0,-0.01, 400.0, 0.2, -150.0, 450.0, -0.01, -150.0, 95.0 };
  vector<double> hi={32, 1.1, 1700.0, 1.8, 150.0, 650.0, 1.1, 250.0, 110.0 };
  
  TCut trigger,goodfit,pid,rmom;
  TTree* ta;
  if(charge<0){
    trigger = _eminustrig;
    goodfit = _eminustq;
    pid = _eminuspid;
    rmom = _eminusrmom;
    ta = _tn;
  } else {
    trigger = _eplustrig;
    goodfit = _eplustq;
    pid = _epluspid;
    rmom = _eplusrmom;
    ta = _tp;
  }

  vector<TCut> cuts(names.size());
  cuts[0] = goodfit+_livegate+_rpitch+_opa+_CRV+pid+rmom;
  cuts[1] = trigger+_livegate+_rpitch+_opa+_CRV+pid+rmom;
  cuts[2] = trigger+goodfit+_rpitch+_opa+_CRV+pid+rmom;
  cuts[3] = trigger+goodfit+_livegate+_opa+_CRV+pid+rmom;
  cuts[4] = trigger+goodfit+_livegate+_rpitch+_CRV+pid+rmom;
  cuts[5] = trigger+goodfit+_livegate+_rpitch+_CRV+pid+rmom;
  cuts[6] = trigger+goodfit+_livegate+_rpitch+_opa+_CRV+rmom;
  cuts[7] = trigger+goodfit+_livegate+_rpitch+_opa+pid+rmom+"bestcrv>=0";
  cuts[8] = trigger+goodfit+_livegate+_rpitch+_opa+_CRV+pid;

  vector<TH1F*> plots;
  for(size_t iplot=0;iplot<names.size();iplot++){
    plots.push_back(new TH1F(names[iplot].c_str(),titles[iplot].c_str(),100,low[iplot],hi[iplot]));
    ta->Project(names[iplot].c_str(),vars[iplot].c_str(),cuts[iplot]);
  }
  _spcan = new TCanvas("spcan","spcan",1000,1000);
  _spcan->Divide(3,3);
  for(size_t iplot=0;iplot<names.size();iplot++){
    _spcan->cd(iplot+1);
    plots[iplot]->Draw();
  }  
}

void TrkAnaPlots::Acc(int ngen,int charge) {
  vector<string> binnames={"All","Trigger", "KF Track fit", "Fit Quality", "Livegate", "Reco pitch", "OPA Rejection", "Upstream", "PID", "CRV Rejection", "Momentum window"};

  unsigned nbins = binnames.size();
  double bmax = nbins-0.5;

  TH1F* acc = new TH1F("acc","Acceptance #times Efficiency;;Cummulative a#times#epsilon",nbins,-0.5,bmax);
  TH1F* racc = new TH1F("racc","Acceptance #times Efficiency;;Relative a#times#epsilon",nbins,-0.5,bmax);
  TH1F* norm = new TH1F("norm","Normalization",nbins,-0.5,bmax);
  TH1F* eff = new TH1F("eff","Cut Efficiency;;#epsilon after all other cuts",nbins,-0.5,bmax);
  TH1F* rej = new TH1F("rej","Cut Rejection;;Fraction left after all other cuts",nbins,-0.5,bmax);
  vector<string> binvals;
  char bval[20];
  for(size_t ibin=0;ibin<nbins;ibin++){
    acc->GetXaxis()->SetBinLabel(ibin+1,binnames[ibin].c_str());
    racc->GetXaxis()->SetBinLabel(ibin+1,binnames[ibin].c_str());
    eff->GetXaxis()->SetBinLabel(ibin+1,binnames[ibin].c_str());
    rej->GetXaxis()->SetBinLabel(ibin+1,binnames[ibin].c_str());
    snprintf(bval,20,"%i",(int)ibin);
    binvals.push_back(string(bval));
  }
  acc->SetStats(0);
  racc->SetStats(0);
  eff->SetStats(0);
  rej->SetStats(0);

  unsigned ibin = 0;
  TCut goodmc("demc.prel>=0");
  TCut trigger,goodfit,pid,rmom;
  TTree* ta;
  if(charge<0){
    trigger = _eminustrig;
    goodfit = _eminustq;
    pid = _eminuspid;
    rmom = _eminusrmom;
    ta = _tn;
  } else {
    trigger = _eplustrig;
    goodfit = _eplustq;
    pid = _epluspid;
    rmom = _eplusrmom;
    ta = _tp;
  }

  ta->Project("acc",binvals[ibin++].c_str(),_pbi*goodmc);
  ta->Project("+acc",binvals[ibin++].c_str(),_pbi*(goodmc+trigger));
  ta->Project("+acc",binvals[ibin++].c_str(),_pbi*(goodmc+trigger+_reco));
  ta->Project("+acc",binvals[ibin++].c_str(),_pbi*(goodmc+trigger+_reco+goodfit));
  ta->Project("+acc",binvals[ibin++].c_str(),_pbi*(goodmc+trigger+_reco+goodfit+_livegate));
  ta->Project("+acc",binvals[ibin++].c_str(),_pbi*(goodmc+trigger+_reco+goodfit+_livegate+_rpitch));
  ta->Project("+acc",binvals[ibin++].c_str(),_pbi*(goodmc+trigger+_reco+goodfit+_livegate+_rpitch+_opa));
  ta->Project("+acc",binvals[ibin++].c_str(),_pbi*(goodmc+trigger+_reco+goodfit+_livegate+_rpitch+_opa+_upstream));
  ta->Project("+acc",binvals[ibin++].c_str(),_pbi*(goodmc+trigger+_reco+goodfit+_livegate+_rpitch+_opa+_upstream+pid));
  ta->Project("+acc",binvals[ibin++].c_str(),_pbi*(goodmc+trigger+_reco+goodfit+_livegate+_rpitch+_opa+_upstream+_CRV+pid));
  ta->Project("+acc",binvals[ibin++].c_str(),_pbi*(goodmc+trigger+_reco+goodfit+_livegate+_rpitch+_opa+_upstream+_CRV+pid+rmom));

  ta->Project("norm",binvals[0].c_str(),_pbi*goodmc);
  ibin = 0;
  ta->Project("+eff",binvals[ibin++].c_str(),(goodmc+trigger+_reco+goodfit+_livegate+_rpitch+_opa+_upstream+_CRV+pid+rmom));
  ta->Project("+eff",binvals[ibin++].c_str(),(goodmc+_reco+goodfit+_livegate+_rpitch+_opa+_upstream+_CRV+pid+rmom));
  ta->Project("+eff",binvals[ibin++].c_str(),(goodmc+trigger+goodfit+_livegate+_rpitch+_opa+_upstream+_CRV+pid+rmom));
  ta->Project("+eff",binvals[ibin++].c_str(),(goodmc+trigger+_reco+_livegate+_rpitch+_opa+_upstream+_CRV+pid+rmom));
  ta->Project("+eff",binvals[ibin++].c_str(),(goodmc+trigger+_reco+goodfit+_rpitch+_opa+_upstream+_CRV+pid+rmom));
  ta->Project("+eff",binvals[ibin++].c_str(),(goodmc+trigger+_reco+goodfit+_livegate+_opa+_upstream+_CRV+pid+rmom));
  ta->Project("+eff",binvals[ibin++].c_str(),(goodmc+trigger+_reco+goodfit+_livegate+_rpitch+_upstream+_CRV+pid+rmom));
  ta->Project("+eff",binvals[ibin++].c_str(),(goodmc+trigger+_reco+goodfit+_livegate+_rpitch+_opa+_CRV+pid+rmom));
  ta->Project("+eff",binvals[ibin++].c_str(),(goodmc+trigger+_reco+goodfit+_livegate+_rpitch+_opa+_upstream+_CRV+rmom));
  ta->Project("+eff",binvals[ibin++].c_str(),(goodmc+trigger+_reco+goodfit+_livegate+_rpitch+_opa+_upstream+pid+rmom));
  ta->Project("+eff",binvals[ibin++].c_str(),(goodmc+trigger+_reco+goodfit+_livegate+_rpitch+_opa+_upstream+_CRV+pid));

  double all = acc->GetBinContent(1);
  if(ngen < 0)
    ngen = all;
  double prev = ngen;
  for(ibin=1;ibin<=nbins;ibin++){
    if(prev > 0.0)
      racc->SetBinContent(ibin,acc->GetBinContent(ibin)/prev);
    else
      racc->SetBinContent(ibin,0.0);
    prev = acc->GetBinContent(ibin);
  }
  racc->SetMaximum(1.1);
  racc->SetMinimum(-0.05);
  acc->Scale(1.0/(float)ngen);
  acc->SetMaximum(1.1);
  acc->SetMinimum(-0.05);
  acc->SetStats(0);
  racc->SetStats(0);
  acc->GetXaxis()->SetLabelSize(0.06);
  racc->GetXaxis()->SetLabelSize(0.06);
  acc->SetMarkerSize(2.0);
  racc->SetMarkerSize(2.0);
  acc->GetYaxis()->SetTitleSize(0.05);
  racc->GetYaxis()->SetTitleSize(0.05);
  eff->GetXaxis()->SetLabelSize(0.06);
  rej->GetXaxis()->SetLabelSize(0.06);
  eff->SetMarkerSize(2.0);
  rej->SetMarkerSize(2.0);
  eff->GetYaxis()->SetTitleSize(0.05);
  rej->GetYaxis()->SetTitleSize(0.05);

  gStyle->SetPaintTextFormat("5.4f");
  _acan = new TCanvas("acan","Acceptance",1200,800);
  _acan->Clear();
  _acan->Divide(1,2);
  _acan->cd(1);
  TPad* tp = (TPad*)_acan->cd(1);
  tp->SetBottomMargin(0.15);
  acc->Draw("histtext0");
  _acan->cd(2);
  tp = (TPad*)_acan->cd(2);
  tp->SetBottomMargin(0.15);
  racc->Draw("histtext0");

  double normval = norm->GetBinContent(1);
  double allval = eff->GetBinContent(1);
  cout << "Found " << normval << "Entries, " << allval << " survive all cuts." << endl;
  for(ibin=1;ibin<=nbins;ibin++){
//    cout << "bin " << ibin << eff->GetXaxis()->GetBinLabel(ibin) <<  " contents = " << eff->GetBinContent(ibin) << endl;
    rej->SetBinContent(ibin,eff->GetBinContent(ibin)/normval);
    if(eff->GetBinContent(ibin)>0.0)
      eff->SetBinContent(ibin,allval/eff->GetBinContent(ibin));
    else
      eff->SetBinContent(ibin,0.0);
  }
  gStyle->SetPaintTextFormat("5.5f");
  _ecan = new TCanvas("ecan","CutEfficiency",1200,1200);
  _ecan->Divide(1,2);
  _ecan->cd(1);
  tp = (TPad*)_ecan->cd(1);
  tp->SetBottomMargin(0.15);
  eff->Draw("histtext0");
  tp = (TPad*)_ecan->cd(2);
  tp->SetBottomMargin(0.15);
  tp->SetLogy();
  rej->Draw("histtext0");
}

void TrkAnaPlots::hitres() {
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
  _tn->Project("hresida","(detsh._doca-detsh._rdrift*detsh._ambig)",_reco+_reco+"detsh._active&&detsh._ambig!=0");
  _tn->Project("hresidna","(detsh._doca-detsh._rdrift*detsh._ambig)",_reco+"detsh._active&&detsh._ambig==0");
  _tn->Project("hresidall","(detsh._doca-detsh._rdrift*detsh._ambig)",_reco+"detsh._active");
  
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
  _tn->Project("hresa","detsh._rdrift-detshmc._dist",_reco+"detsh._active&&detsh._ambig!=0");
  _tn->Project("hresna","detsh._rdrift-detshmc._dist",_reco+"detsh._active&&detsh._ambig==0");
  _tn->Project("hresall","detsh._rdrift-detshmc._dist",_reco+"detsh._active");
  _rescan = new TCanvas("rescan","rescan",800,800);
  _rescan->Divide(1,2);
  _rescan->cd(1);
  hresidst->Draw("h");
  hresidall->Fit("gaus","","sames");
  TLegend* rleg = new TLegend(0.15,0.6,0.4,0.85);
  rleg->AddEntry(hresida,"Resolved Ambiguity","f");
  rleg->AddEntry(hresidna,"Null Ambiguity","f");
  rleg->Draw();
  _rescan->cd(2);
  hresst->Draw("h");
  hresall->Fit("gaus","","sames");
}

void TrkAnaPlots::wpull() {
  TH1F* swp = new TH1F("swp","Final Fit Wire Position Pull",100,-25,25);
  TH1F* uwp = new TH1F("uwp","Final Fit Wire Position Pull",100,-25,25);
  TH1F* rwp = new TH1F("rwp","Final Fit Wire Position Pull",100,-25,25);
  swp->SetLineColor(kGreen);
  uwp->SetLineColor(kRed);
  rwp->SetLineColor(kBlue);
  _tn->Project("swp","(detsh._wdist-detsh._hlen)/detsh._werr",_reco+"detsh._active&&detshmc._rel._rel==0");
  _tn->Project("uwp","(detsh._wdist-detsh._hlen)/detsh._werr",_reco+"detsh._active&&detshmc._rel._rel<0");
  _tn->Project("rwp","(detsh._wdist-detsh._hlen)/detsh._werr",_reco+"detsh._active&&detshmc._rel._rel>0");
  _wpcan = new TCanvas("wpcan","wpcan",800,800);
  _wpcan->Divide(2,2);
  _wpcan->cd(1);
  gPad->SetLogy();
  swp->Fit("gaus");
  _wpcan->cd(2);
  uwp->Draw();
  _wpcan->cd(3);
  gPad->SetLogy();
  rwp->Fit("gaus");
  TLegend* wleg = new TLegend(0.5,0.5,0.9,0.9);
  wleg->AddEntry(swp,"Primary Hit","L");
  wleg->AddEntry(uwp,"Unrelated Hit","L");
  wleg->AddEntry(rwp,"Related Hit","L");
  _wpcan->cd(4);
  wleg->Draw();

}

void TrkAnaPlots::Ambig() {
  gStyle->SetOptStat(1111);

  TCut ghit("detshmc._rel._rel==0");
  TCut delta("detshmc._rel._rel>0");
  TCut bkg("detshmc._rel._rel<0");
  TCut gambig("detshmc._ambig==detsh._ambig");
  TCut bambig("detshmc._ambig!=detsh._ambig&&detsh._ambig!=0");
  TCut nambig("detsh._ambig==0");
  TCut active("detsh._active>0");
// apply requested cuts

  TCut goodtrk(_reco+"de.nactive>20");


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
  rdg->SetMaximum(1.1);
  rdg->SetMinimum(-0.1);
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

  _tn->Project("rdg","detshmc._dist",goodtrk+active+gambig+ghit);
  _tn->Project("rdn","detshmc._dist",goodtrk+active+nambig+ghit);
  _tn->Project("rdb","detshmc._dist",goodtrk+active+bambig+ghit);
  _tn->Project("rda","detshmc._dist",goodtrk+active);
  _tn->Project("rdi","detshmc._dist",goodtrk+ghit+(!active));
  _tn->Project("rdd","detshmc._dist",goodtrk+active+delta);
  _tn->Project("rdf","detshmc._dist",goodtrk+active+bkg);
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
  _tn->Project("momres0","deent.mom-sqrt(demcent.momx^2+demcent.momy^2+demcent.momz^2)",goodtrk);
  _tn->Project("momres1","deent.mom-sqrt(demcent.momx^2+demcent.momy^2+demcent.momz^2)",goodtrk&&"de.status==1");
  _tn->Project("momres2","deent.mom-sqrt(demcent.momx^2+demcent.momy^2+demcent.momz^2)",goodtrk&&"de.status==2");

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
  _tn->Project("afg","deent.mom-sqrt(demcent.momx^2+demcent.momy^2+demcent.momz^2)",goodtrk+active+gambig);
  _tn->Project("afn","deent.mom-sqrt(demcent.momx^2+demcent.momy^2+demcent.momz^2)",goodtrk+active+nambig);
  _tn->Project("afb","deent.mom-sqrt(demcent.momx^2+demcent.momy^2+demcent.momz^2)",goodtrk+active+bambig);
  _tn->Project("afa","deent.mom-sqrt(demcent.momx^2+demcent.momy^2+demcent.momz^2)",goodtrk+active);
  afg->Divide(afa);
  afn->Divide(afa);
  afb->Divide(afa);
  afg->SetMinimum(0.0);
  afg->SetMaximum(1.1);

  _ambigcan = new TCanvas("ambigcan","Hit Ambiguity",1200,800);
  _ambigcan->Divide(2,2);


  _ambigcan->cd(1);
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

  _ambigcan->cd(2);

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

  _ambigcan->cd(3);
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


  _ambigcan->cd(4);
  afg->Draw();
  afn->Draw("same");
  afb->Draw("same");
  TLegend* fleg = new TLegend(0.16,0.35,0.625,0.6);
  fleg->AddEntry(afg,"Correct Ambiguity","l");
  fleg->AddEntry(afn,"No Ambiguity","l");
  fleg->AddEntry(afb,"Incorrect Ambiguity","l");
  fleg->Draw();

  _ambigcan->cd(0);
}

void TrkAnaPlots::Resid() {

  TCut delta("detshmc._proc==17");
  TCut primary("detshmc._gen==2");
  TCut gambig("detshmc._ambig==detsh._ambig&&detsh._ambig!=0");
  TCut bambig("(detshmc._ambig!=detsh._ambig)&&detsh._ambig!=0");
  TCut nambig("detsh._ambig==0");
  TCut active("de.status==1 && detsh._active>0");
  TCut reco(_reco);
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

  _tn->Project("rdg","detshmc._dist",reco+mcsel+active+gambig+primary);
  _tn->Project("rdb","detshmc._dist",reco+mcsel+active+bambig+primary);
  _tn->Project("rdn","detshmc._dist",reco+mcsel+active+nambig+primary);

  _tn->Project("rpullg","(detsh._doca-detsh._rdrift*detsh._ambig)/detsh._residerr",reco+mcsel+active+gambig+primary);
  _tn->Project("rpullb","(detsh._doca-detsh._rdrift*detsh._ambig)/detsh._residerr",reco+mcsel+active+bambig+primary);
  _tn->Project("rpulln","(detsh._doca-detsh._rdrift*detsh._ambig)/detsh._residerr",reco+mcsel+active+nambig+primary);
  _tn->Project("rpulld","(detsh._doca-detsh._rdrift*detsh._ambig)/detsh._residerr",reco+mcsel+active+delta);

  _residcan = new TCanvas("residcan","Residuals",1200,800);
  _residcan->Divide(2,1);
  _residcan->cd(1);
  rdg->Draw();
  rdb->Draw("same");
  rdn->Draw("same");

  TLegend* leg = new TLegend(0.3,0.3,0.8,0.5);
  leg->AddEntry(rpullg,"Correct ambiguity","l");
  leg->AddEntry(rpullb,"Incorrect ambiguity","l");
  leg->AddEntry(rpulln,"No ambiguity assigned","l");
  leg->Draw();

  TPad* ppad = dynamic_cast<TPad*>(_residcan->cd(2));
  ppad->Divide(1,4);
  ppad->cd(1);
  rpullg->Fit("gaus");
  ppad->cd(2);
  rpullb->Fit("gaus");
  ppad->cd(3);
  rpulln->Fit("gaus");
  ppad->cd(4);
  rpulld->Fit("gaus");

  _residcan->cd(0);

}

void TrkAnaPlots::Con() {
  TCut mcsel("sqrt(demcent.momx^2+demcent.momy^2+demcent.momz^2)>100.0");

  TH1F* con1 = new TH1F("con1","#chi^{2} fit consistency",500,0.0,1.0);
  TH1F* con2 = new TH1F("con2","#chi^{2} fit consistency",500,0.0,1.0);
  TH1F* lcon1 = new TH1F("lcon1","log_{10}(#chi^{2}) fit consistency",100,-10,0);
  TH1F* lcon2 = new TH1F("lcon2","log_{10}(#chi^{2}) fit consistency",100,-10,0);
  con1->SetLineColor(kBlue);
  con2->SetLineColor(kRed);
  lcon1->SetLineColor(kBlue);
  lcon2->SetLineColor(kRed);

  _tn->Project("con1","de.fitcon",mcsel+"de.status==1");
  _tn->Project("con2","de.fitcon",mcsel+"de.status==2");
  _tn->Project("lcon1","log10(de.fitcon)",mcsel+"de.status==1");
  _tn->Project("lcon2","log10(de.fitcon)",mcsel+"de.status==2");

  _fcan = new TCanvas("fcan","fit consistency",500,800);
  _fcan->Clear();
  _fcan->Divide(1,2);
  _fcan->cd(1);
  con1->Draw();
  con2->Draw("same");
  _fcan->cd(2);
  lcon1->Draw();
  lcon2->Draw("same");

  TLegend* leg = new TLegend(0.1,0.6,0.4,0.8);
  leg->AddEntry(con1,"Fully Converged Fit","l");
  leg->AddEntry(con2,"Unconverged Fit","l");
  leg->Draw();

}

void TrkAnaPlots::TrkQual(const char* extra) {
  TCut ecut(extra);
  TH1F* tq = new TH1F("tq","TrkQual;TrkQual MVA Output",203,-0.01,1.01);
  TH1F* tqtch = new TH1F("tqtch","TrkQual;TrkQual MVA Output",203,-0.01,1.01);
  TH1F* tqntch = new TH1F("tqntch","TrkQual;TrkQual MVA Output",203,-0.01,1.01);
  tq->SetLineColor(kBlack);
  tqtch->SetLineColor(kBlue);
  tqntch->SetLineColor(kRed);
  tq->SetStats(0);
  tqtch->SetStats(0);
  tqntch->SetStats(0);
  _tn->Project("tq","dequal.TrkQualDeM",ecut);
  _tn->Project("tqtch","dequal.TrkQualDeM","detch.active"+ecut);
  _tn->Project("tqntch","dequal.TrkQualDeM","!detch.active"+ecut);

  double* tqintarray = tq->GetIntegral();

  TH1F* tqint = new TH1F("tqint","TrkQual Cut Efficiency;TrkQual Cut;Efficiency",tq->GetNbinsX(),tq->GetXaxis()->GetXmin(), tq->GetXaxis()->GetXmax());
  tqint->SetStats(0);
  tqint->SetFillColor(kBlack);
  for(size_t ibin=0;ibin < (size_t)tq->GetNbinsX();ibin++)
    tqint->SetBinContent(ibin+1, 1.0-tqintarray[ibin]);
  _tqcan = new TCanvas("tqcan","TrkQual",600,600);
  _tqcan->Divide(1,2);
  _tqcan->cd(1);
  gPad->SetLogy();
  tq->Draw();
  tqtch->Draw("same");
  tqntch->Draw("same");
  TLegend* tqleg = new TLegend(0.5,0.7,0.8,0.9);
  tqleg->AddEntry(tq,"All","L");
  tqleg->AddEntry(tqtch,"TrkCaloHit","L");
  tqleg->AddEntry(tqntch,"No TrkCaloHit","L");
  tqleg->Draw();
  _tqcan->cd(2);
  tqint->Draw();
}

void TrkAnaPlots::TrkQualRes(float tqcut) {
  TF1* dscb = new TF1("dscb",fnc_dscb,-2.0,4,7);
  dscb->SetParName(0,"Norm");
  dscb->SetParName(1,"x0");
  dscb->SetParName(2,"sigma");
  dscb->SetParName(3,"ANeg");
  dscb->SetParName(4,"PNeg");
  dscb->SetParName(5,"APos");
  dscb->SetParName(6,"PPos");
  dscb->SetNpx(1000);
  dscb->SetParLimits(3,0.0,50.0);
  dscb->SetParLimits(4,1.0,50.0);
  dscb->SetParLimits(5,0.0,50.0);
  dscb->SetParLimits(6,1.0,50.0);

   TH1F* goodf = new TH1F("goodf","Momentum Resolution;Reco - True Momentum (MeV/c)",100,-4,4);
  TH1F* badf = new TH1F("badf","Momentum Resolution;Reco - True Momentum (MeV/c)",100,-4,4);
  goodf->SetLineColor(kBlue);
  goodf->SetMarkerColor(kBlue);
  goodf->SetLineWidth(2);
  goodf->SetFillColor(kBlue);
  goodf->SetFillStyle(3004);
  badf->SetLineColor(kRed);
  badf->SetFillColor(kRed);
  badf->SetMarkerColor(kRed);
  badf->SetLineWidth(2);
  badf->SetFillStyle(3005);
  goodf->SetStats(0);
  badf->SetStats(0);
  goodf->Sumw2();
  badf->Sumw2();
  char tqcutgc[40];
  char tqcutbc[40];
  snprintf(tqcutgc,40,"dequal.TrkQualDeM>%f",tqcut);
  snprintf(tqcutbc,40,"dequal.TrkQualDeM>0");
  TCut tqcutg(tqcutgc);
  TCut tqcutb(tqcutbc);
  TCut reco(_reco);
  TCut mcsel("sqrt(demcent.momx^2+demcent.momy^2+demcent.momz^2)>100.0");
  _tn->Project("goodf","deent.mom-sqrt(demcent.momx^2+demcent.momy^2+demcent.momz^2)",(reco+mcsel+tqcutg)*"_pbi.PBIWeight");
  _tn->Project("badf","deent.mom-sqrt(demcent.momx^2+demcent.momy^2+demcent.momz^2)",(reco+mcsel)*"_pbi.PBIWeight");
  
  _tqrcan = new TCanvas("tqrcan","TrkQualRes",1000,800);
  TLegend* leg = new TLegend(0.6,0.6,0.9,0.9);
  char ltitle[40];
  gPad->SetLogy();

  double integral = badf->GetEntries()*badf->GetBinWidth(1);
  dscb->SetParameters(integral,0.0,0.15,1.0,4.5,1.2,10.0);
  dscb->SetLineColor(kRed);
  badf->Fit("dscb","LIR");
  double sigma = dscb->GetParameter(2);
  double sigerr = dscb->GetParError(2);
  double ppos = dscb->GetParameter(6);
  double pposerr = dscb->GetParError(6);
  leg->AddEntry(badf,"No TrkQual Cut","L");
  snprintf(ltitle,40,"#sigma = %3.1f #pm %2.2f KeV/c",sigma*1000,sigerr*1000);
  leg->AddEntry(badf,ltitle,"L");
  snprintf(ltitle,40,"PPos = %3.1f #pm %2.2f",ppos,pposerr);
  leg->AddEntry(badf,ltitle,"L");

  integral = goodf->GetEntries()*badf->GetBinWidth(1);
  dscb->SetParameters(integral,0.0,0.15,1.0,4.5,1.2,10.0);
  dscb->SetLineColor(kBlue);
  goodf->Fit("dscb","LIR","sameLF2");
  snprintf(ltitle,40,"TrkQual>%3.2f, f = %5.3f",tqcut,goodf->GetEntries()/float(badf->GetEntries()));
  leg->AddEntry(goodf,ltitle,"L");
  sigma = dscb->GetParameter(2);
  sigerr = dscb->GetParError(2);
  ppos = dscb->GetParameter(6);
  pposerr = dscb->GetParError(6);
  snprintf(ltitle,40,"#sigma = %3.1f #pm %2.2f KeV/c",sigma*1000,sigerr*1000);
  leg->AddEntry(goodf,ltitle,"L");
  snprintf(ltitle,40,"PPos = %3.1f #pm %2.2f",ppos,pposerr);
  leg->AddEntry(goodf,ltitle,"L");

  leg->Draw();
}

void TrkAnaPlots::StrawMat() {
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

  _tn->Project("addmatfrac","(de.nmatactive-de.nactive)/de.nmatactive",_reco);
  _tn->Project("matvshit","de.nmatactive:de.nactive",_reco);
  _tn->Project("naddmat","de.nmatactive-de.nactive",_reco);
  _tn->Project("adddoca","detsm._doca",_reco+"detsm._active&&(!detsm._thita)");
  _tn->Project("hitdoca","detsm._doca",_reco+"detsm._thita");
  _tn->Project("addstraw","detsm._straw",_reco+"detsm._active&&(!detsm._thita)");
  _tn->Project("hitstraw","detsm._straw",_reco+"detsm._thita");

  _tn->Project("lofracres","deent.mom-sqrt(demcent.momx^2+demcent.momy^2+demcent.momz^2)",_reco+"(de.nmatactive-de.nactive)/de.nmatactive<0.1");
  _tn->Project("hifracres","deent.mom-sqrt(demcent.momx^2+demcent.momy^2+demcent.momz^2)",_reco+"(de.nmatactive-de.nactive)/de.nmatactive>0.1");

  TLegend* leg = new TLegend(0.6,0.7,0.9,.9);
  leg->AddEntry(hitdoca,"Hit Straw","L");
  leg->AddEntry(adddoca,"Added Straw","L");

  _mcan = new TCanvas("mcan","mcan",1200,800);
  _mcan->Divide(3,2);
  _mcan->cd(1);
  matvshit->Draw("colorz");
  _mcan->cd(2);
  naddmat->Draw();
  _mcan->cd(3);
  addmatfrac->Draw();
  _mcan->cd(4);
  hitdoca->Draw();
  adddoca->Draw("same");
  leg->Draw();
  _mcan->cd(5);
  hitstraw->Draw();
  addstraw->Draw("same");
  leg->Draw();
  _mcan->cd(6);
  lofracres->Draw();
  hifracres->Draw("same");
  TLegend* mleg = new TLegend(0.6,0.7,0.9,0.9);
  mleg->AddEntry(lofracres,"N_{added}/N<0.1","L");
  mleg->AddEntry(hifracres,"N_{added}/N>0.1","L");
  mleg->Draw();
}

void TrkAnaPlots::TrkCaloHit(float tqcut,int pdg) {
  char cstring[100];
  snprintf(cstring,100,"detch.active&&dequal.TrkQualDeM>%f&&abs(demc.pdg)==%i&&demcxit.momz>0",tqcut,pdg);
  TCut goodtrkcalo(cstring);
  snprintf(cstring,100,"(!detch.active)&&dequal.TrkQualDeM>%f&&abs(demc.pdg)==%i&&demcxit.momz>0",tqcut,pdg);
  TCut goodtrk(cstring);
  TCut disk0("detch.disk==0");
  TCut disk1("detch.disk==1");
  TH1F* clen0 = new TH1F("clen0","TrkCaloHit POCA Crystal Depth;Depth(mm)",200,-50,250);
  TH1F* clen1 = new TH1F("clen1","TrkCaloHit POCA Crystal Depth;Depth(mm)",200,-50,250);
  TH1F* cdoca0 = new TH1F("cdoca0","TrkCaloHit DOCA;DOCA (mm)",200,-150,150);
  TH1F* cdoca1 = new TH1F("cdoca1","TrkCaloHit DOCA;DOCA (mm)",200,-150,150);
  TH1F* cdt0 = new TH1F("cdt0","TrkCaloHit #Delta t;TCH t_{0} - CaloCluster Time",100,-5,5);
  TH1F* cdt1 = new TH1F("cdt1","TrkCaloHit #Delta t;TCH t_{0} - CaloCluster Time",100,-5,5);
  TH1F* ep0 = new TH1F("ep0","TrkCaloHit E/P",100,0.0,1.25);
  TH1F* ep1 = new TH1F("ep1","TrkCaloHit E/P",100,0.0,1.25);
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
 
  _tn->Project("clen0","detch.clen",goodtrkcalo&&disk0);
  _tn->Project("clen1","detch.clen",goodtrkcalo&&disk1);
  _tn->Project("cdoca0","detch.doca",goodtrkcalo&&disk0);
  _tn->Project("cdoca1","detch.doca",goodtrkcalo&&disk1);
  _tn->Project("cdt0","detch.t0-detch.ctime",goodtrkcalo&&disk0);
  _tn->Project("cdt1","detch.t0-detch.ctime",goodtrkcalo&&disk1);
  _tn->Project("ep0","detch.edep/sqrt(detch.momx^2+detch.momy^2+detch.momz^2)",goodtrkcalo&&disk0);
  _tn->Project("ep1","detch.edep/sqrt(detch.momx^2+detch.momy^2+detch.momz^2)",goodtrkcalo&&disk1);
  _tn->Project("tdir0","(detch.POCAx*detch.momx+detch.POCAy*detch.momy)/sqrt(detch.POCAx^2+detch.POCAy^2)/sqrt(detch.momx^2+detch.momy^2+detch.momz^2)",goodtrkcalo&&disk0);
  _tn->Project("tdir1","(detch.POCAx*detch.momx+detch.POCAy*detch.momy)/sqrt(detch.POCAx^2+detch.POCAy^2)/sqrt(detch.momx^2+detch.momy^2+detch.momz^2)",goodtrkcalo&&disk1);
  _tn->Project("pr0","sqrt(detch.POCAx^2+detch.POCAy^2)",goodtrkcalo&&disk0);
  _tn->Project("pr1","sqrt(detch.POCAx^2+detch.POCAy^2)",goodtrkcalo&&disk1);

  
  TLegend* tchleg = new TLegend(0.6,0.7,0.9,0.9);
  tchleg->AddEntry(clen0,"Disk 0","L");
  tchleg->AddEntry(clen1,"Disk 1","L");

  _tchcan = new TCanvas("tchcan","TrkCaloHit",800,600);
  _tchcan->Divide(3,2);
  _tchcan->cd(1);
  clen0->Draw();
  clen1->Draw("same");
  tchleg->Draw();
  _tchcan->cd(2);
  cdoca0->Draw();
  cdoca1->Draw("same");
  _tchcan->cd(3);
  cdt0->Draw();
  cdt1->Draw("same");
  _tchcan->cd(4);
  ep0->Draw();
  ep1->Draw("same");
  _tchcan->cd(5);
  tdir0->Draw();
  tdir1->Draw("same");
  _tchcan->cd(6);
  pr0->Draw();
  pr1->Draw("same");
  

  TCut goodclen("detch.clen>0&&detch.clen<150.0");
  TCut badclen("detch.clen>150.0&&detch.clen<250.0");
  TH2F* rad0 = new TH2F("rad0","POCA Radius vs Crystal Depth, Disk 0;Depth (mm);Radius (mm)",100,-50,250,100,360,650);
  TH2F* rad1 = new TH2F("rad1","POCA Radius vs Crystal Depth, Disk 1;Depth (mm);Radius (mm)",100,-50,250,100,360,650);
  TH2F* dot0 = new TH2F("dot0","Track Direction vs Crystal Depth, Disk 0;Depth (mm);#hat{t}#bullet#hat{#rho}",100,-50,250,100,-1,1);
  TH2F* dot1 = new TH2F("dot1","Track Direction vs Crystal Depth, Disk 1;Depth (mm);#hat{t}#bullet#hat{#rho}",100,-50,250,100,-1,1);
  TH2F* dotr0 = new TH2F("dotr0","Track Direction vs POCA Radius, Disk 0;Radius (mm);#hat{t}#bullet#hat{#rho}",100,360,650,100,-1,1);
  TH2F* dotr1 = new TH2F("dotr1","Track Direction vs POCA Radius, Disk 1;Radius (mm);#hat{t}#bullet#hat{#rho}",100,360,650,100,-1,1);
  TH2F* doca0 = new TH2F("doca0","DOCA vs Crystal Depth, Disk 0;Depth (mm);DOCA (mm)",100,-50,250,100,-100,100);
  TH2F* doca1 = new TH2F("doca1","DOCA vs Crystal Depth, Disk 1;Depth (mm);DOCA(mm)",100,-50,250,100,-100,100);
  TH2F* eopd0 = new TH2F("eopd0","E/P vs Crystal Depth, Disk 0;Depth (mm);E/P",100,-50,250,100,0.0,1.25);
  TH2F* eopd1 = new TH2F("eopd1","E/P vs Crystal Depth, Disk 1;Depth (mm);E/P",100,-50,250,100,0.0,1.25);
  TH2F* eopdir0 = new TH2F("eopdir0","E/P vs Track Direction, Disk 0;#hat{t}#bullet#hat{#rho};E/P",100,-1,1,100,0.0,1.25);
  TH2F* eopdir1 = new TH2F("eopdir1","E/P vs Track Direction, Disk 1;#hat{t}#bullet#hat{#rho};E/P",100,-1,1,100,0.0,1.25);
  TH2F* eopr0 = new TH2F("eopr0","E/P vs POCA Radius, Disk 0;Radius (mm);E/P",100,360,650,100,0.0,1.25);
  TH2F* eopr1 = new TH2F("eopr1","E/P vs POCA Radius, Disk 1;Radius (mm);E/P",100,360,650,100,0.0,1.25);

  TProfile* peopd0 = new TProfile("peopd0","E/P vs Crystal Depth, Disk 0;Depth (mm);E/P",100,-50,250,0.0,1.25);
  TProfile* peopd1 = new TProfile("peopd1","E/P vs Crystal Depth, Disk 1;Depth (mm);E/P",100,-50,250,0.0,1.25);
 
  rad0->SetStats(0);
  rad1->SetStats(0);
  dot0->SetStats(0);
  dot1->SetStats(0);
  doca0->SetStats(0);
  doca1->SetStats(0);
  eopd0->SetStats(0);
  eopd1->SetStats(0);
  eopdir0->SetStats(0);
  eopdir1->SetStats(0);
  eopr0->SetStats(0);
  eopr1->SetStats(0);
  _tn->Project("rad0","sqrt(detch.POCAx^2+detch.POCAy^2):detch.clen",goodtrkcalo&&disk0);
  _tn->Project("rad1","sqrt(detch.POCAx^2+detch.POCAy^2):detch.clen",goodtrkcalo&&disk1);
  _tn->Project("dot0","(detch.POCAx*detch.momx+detch.POCAy*detch.momy)/sqrt(detch.POCAx^2+detch.POCAy^2)/sqrt(detch.momx^2+detch.momy^2+detch.momz^2):detch.clen",goodtrkcalo&&disk0);
  _tn->Project("dot1","(detch.POCAx*detch.momx+detch.POCAy*detch.momy)/sqrt(detch.POCAx^2+detch.POCAy^2)/sqrt(detch.momx^2+detch.momy^2+detch.momz^2):detch.clen",goodtrkcalo&&disk1);
  _tn->Project("dotr0","(detch.POCAx*detch.momx+detch.POCAy*detch.momy)/sqrt(detch.POCAx^2+detch.POCAy^2)/sqrt(detch.momx^2+detch.momy^2+detch.momz^2):sqrt(detch.POCAx^2+detch.POCAy^2)",goodtrkcalo&&disk0);
  _tn->Project("dotr1","(detch.POCAx*detch.momx+detch.POCAy*detch.momy)/sqrt(detch.POCAx^2+detch.POCAy^2)/sqrt(detch.momx^2+detch.momy^2+detch.momz^2):sqrt(detch.POCAx^2+detch.POCAy^2)",goodtrkcalo&&disk1);
  _tn->Project("eopd0","detch.edep/sqrt(detch.momx^2+detch.momy^2+detch.momz^2):detch.clen",goodtrkcalo&&disk0);
  _tn->Project("eopd1","detch.edep/sqrt(detch.momx^2+detch.momy^2+detch.momz^2):detch.clen",goodtrkcalo&&disk1);
  _tn->Project("eopdir0","detch.edep/sqrt(detch.momx^2+detch.momy^2+detch.momz^2):(detch.POCAx*detch.momx+detch.POCAy*detch.momy)/sqrt(detch.POCAx^2+detch.POCAy^2)/sqrt(detch.momx^2+detch.momy^2+detch.momz^2)",goodtrkcalo&&disk0);
  _tn->Project("eopdir1","detch.edep/sqrt(detch.momx^2+detch.momy^2+detch.momz^2):(detch.POCAx*detch.momx+detch.POCAy*detch.momy)/sqrt(detch.POCAx^2+detch.POCAy^2)/sqrt(detch.momx^2+detch.momy^2+detch.momz^2)",goodtrkcalo&&disk1);
  _tn->Project("eopr0","detch.edep/sqrt(detch.momx^2+detch.momy^2+detch.momz^2):sqrt(detch.POCAx^2+detch.POCAy^2)",goodtrkcalo&&disk0);
  _tn->Project("eopr1","detch.edep/sqrt(detch.momx^2+detch.momy^2+detch.momz^2):sqrt(detch.POCAx^2+detch.POCAy^2)",goodtrkcalo&&disk1);
  _tn->Project("doca0","detch.doca:detch.clen",goodtrkcalo&&disk0);
  _tn->Project("doca1","detch.doca:detch.clen",goodtrkcalo&&disk1);

  _tn->Project("peopd0","detch.edep/sqrt(detch.momx^2+detch.momy^2+detch.momz^2):detch.clen",goodtrkcalo&&disk0);
  _tn->Project("peopd1","detch.edep/sqrt(detch.momx^2+detch.momy^2+detch.momz^2):detch.clen",goodtrkcalo&&disk1);

  _tch0can = new TCanvas("tch0can","TrkCaloHit Disk 0", 1000, 800);
  _tch0can->Divide(3,2);
  _tch0can->cd(1);
  gPad->SetLogz();
  rad0->Draw("colorz");
  _tch0can->cd(2);
  gPad->SetLogz();
  dot0->Draw("colorz");
  _tch0can->cd(3);
  gPad->SetLogz();
  eopd0->Draw("colorz");
  peopd0->Draw("same");
  _tch0can->cd(4);
  gPad->SetLogz();
  doca0->Draw("colorz");
  _tch0can->cd(5);
  gPad->SetLogz();
  eopdir0->Draw("colorz");
  _tch0can->cd(6);
  gPad->SetLogz();
  eopr0->Draw("colorz");

  _tch1can = new TCanvas("tch1can","TrkCaloHit Disk 1", 1000, 800);
  _tch1can->Divide(3,2);
  _tch1can->cd(1);
  gPad->SetLogz();
  rad1->Draw("colorz");
  _tch1can->cd(2);
  gPad->SetLogz();
  dot1->Draw("colorz");
  _tch1can->cd(3);
  gPad->SetLogz();
  eopd1->Draw("colorz");
  peopd1->Draw("same");
  _tch1can->cd(4);
  gPad->SetLogz();
  doca1->Draw("colorz");
  _tch1can->cd(5);
  gPad->SetLogz();
  eopdir1->Draw("colorz");
  _tch1can->cd(6);
  gPad->SetLogz();
  eopr1->Draw("colorz");

  TH2F* dtvsclen = new TH2F("dtvsclen","T_{calo} - T_{trk} vs Cluster Depth;Depth (mm); #Delta t (ns)",50,-100,300,50,-2.0,3.0);
  TProfile* pdtvsclen = new TProfile("pdtvsclen","T_{calo} - T_{trk} vs Cluster Depth;Depth (mm); #Delta t (ns)",50,-100,300,-2.0,3.0);
  dtvsclen->SetStats(0);
  _tn->Project("dtvsclen","detch.ctime-detch.t0:detch.clen",goodtrkcalo);
  _tn->Project("pdtvsclen","detch.ctime-detch.t0:detch.clen",goodtrkcalo);

  _dtchtcan = new TCanvas("dtchtcan","TCH timing",800,800);
  gPad->SetLogz();
  dtvsclen->Draw("colorz");
  pdtvsclen->Fit("pol1","","same");

  TH1F* dt0hiE = new TH1F("dt0hiE","T_{0} Resolution;Reco t_{0} - MC t_{0} (nsec)",100,-5,5);
  TH1F* dt0loE = new TH1F("dt0loE","T_{0} Resolution;Reco t_{0} - MC t_{0} (nsec)",100,-5,5);
  TH1F* dt0nocal = new TH1F("dt0nocal","T_{0} Resolution;Reco t_{0} - MC t_{0} (nsec)",100,-5,5);
  dt0hiE->SetStats(0);
  dt0loE->SetStats(0);
  dt0nocal->SetStats(0);
  dt0hiE->SetLineColor(kBlue);
  dt0loE->SetLineColor(kGreen);
  dt0nocal->SetLineColor(kBlack);
  _tn->Project("dt0hiE","de.t0-demcmid.t0",goodtrkcalo&&"detch.edep>50.0");
  _tn->Project("dt0loE","de.t0-demcmid.t0",goodtrkcalo&&"detch.edep<50.0");
  _tn->Project("dt0nocal","de.t0-demcmid.t0",goodtrk);

  TCanvas* dt0can = new TCanvas("dt0can","t0 resolution",800,800);
  gPad->SetLogy();
  dt0hiE->Fit("gaus");
  dt0loE->Fit("gaus","","same");
  dt0nocal->Fit("gaus","","same");
  TLegend* t0leg = new TLegend(0.6,0.6,0.9,0.9);
  t0leg->AddEntry(dt0hiE,"ECalo > 50 MeV/c","L");
  t0leg->AddEntry(dt0loE,"ECalo < 50 MeV/c","L");
  t0leg->AddEntry(dt0nocal,"No TrkCaloHit","L");
  t0leg->Draw();
}

void TrkAnaPlots::TrkCaloHitMC() {
  TCut goodtrkcalo("dequal.TrkQualDeM>0.6&&detch.active");
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
  _tn->Project("clen0m","detch.clen",goodtrkcalo&&disk0&&tcm);
  _tn->Project("clen0n","detch.clen",goodtrkcalo&&disk0&&tcn);
  _tn->Project("clen1m","detch.clen",goodtrkcalo&&disk1&&tcm);
  _tn->Project("clen1n","detch.clen",goodtrkcalo&&disk1&&tcn);

  _tn->Project("cdoca0m","detch.doca",goodtrkcalo&&disk0&&tcm);
  _tn->Project("cdoca0n","detch.doca",goodtrkcalo&&disk0&&tcn);
  _tn->Project("cdoca1m","detch.doca",goodtrkcalo&&disk1&&tcm);
  _tn->Project("cdoca1n","detch.doca",goodtrkcalo&&disk1&&tcn);

  TLegend* tchmcleg = new TLegend(0.5,0.7,0.9,0.9);
  tchmcleg->AddEntry(clen0m,"Trk-Calo MC Match","L");
  tchmcleg->AddEntry(clen0n,"No MC Match","L");

  _tchmccan = new TCanvas("tchmc","TrkCaloHitMC",800,800);
  _tchmccan->Divide(2,2);
  _tchmccan->cd(1);
  gPad->SetLogy();
  clen0m->Draw();
  clen0n->Draw("same");
  _tchmccan->cd(2);
  gPad->SetLogy();
  clen1m->Draw();
  clen1n->Draw("same");
  tchmcleg->Draw();
  _tchmccan->cd(3);
  gPad->SetLogy();
  cdoca0m->Draw();
  cdoca0n->Draw("same");
  _tchmccan->cd(4);
  gPad->SetLogy();
  cdoca1m->Draw();
  cdoca1n->Draw("same");
}

void TrkAnaPlots::t0() {
  TCut goodtrkcalo("dequal.TrkQualDeM>0.6&&detch.active");
  TCut goodtrknocalo("dequal.TrkQualDeM>0.6&&!detch.active");
  TCut disk0("detch.disk==0");
  TCut disk1("detch.disk==1");
  TH1F* t00 = new TH1F("t00","Track Fit t_{0} Resolution, TrkCaloHit;t_{0} reco - t_{0} MC (ns)",100,-5,5);
  TH1F* t01 = new TH1F("t01","Track Fit t_{0} Resolution, NoTrkCaloHit;t_{0} reco - t_{0} MC (ns)",100,-5,5);
//  t00->SetStats(0);
//  t01->SetStats(0);
//  t00->SetLineColor(kRed);
//  t01->SetLineColor(kBlue);
  _tn->Project("t00","de.t0-fmod(demcmid.t0,1695)",goodtrkcalo);
  _tn->Project("t01","de.t0-fmod(demcmid.t0,1695)",goodtrknocalo);
  TLegend* tchleg = new TLegend(0.6,0.7,0.9,0.9);
  tchleg->AddEntry(t00,"TrkCaloHit","L");
  tchleg->AddEntry(t01,"No TrkCaloHit","L");
  _t0can = new TCanvas("t0can","t0can",600,400);
  _t0can->Divide(2,1);
  _t0can->cd(1);
  t00->Draw();
  _t0can->cd(2);
  t01->Draw();
}

void TrkAnaPlots::Eff(unsigned norm, double plo, double phi, int q) {
  TH1F* allrec = new TH1F("allrec","Reco Fraction vs Generated Momentum",100,plo,phi);
  TH1F* tqrec = new TH1F("tqrec","Reco Fraction vs Generated Momentum",100,plo,phi);
  TH1F* t0rec = new TH1F("t0rec","Reco Fraction vs Generated Momentum",100,plo,phi);
  TH1F* tprtrig = new TH1F("tprtrig","Reco Fraction vs Generated Momentum",100,plo,phi);
  TH1F* cprtrig = new TH1F("cprtrig","Reco Fraction vs Generated Momentum",100,plo,phi);
  TH1F* trktrig = new TH1F("trktrig","Reco Fraction vs Generated Momentum",100,plo,phi);
  TH1F* alltrig = new TH1F("alltrig","Reco Fraction vs Generated Momentum",100,plo,phi);
  TH1F* cctrig = new TH1F("cctrig","Reco Fraction vs Generated Momentum",100,plo,phi);

  TCut t0cut("de.t0>700");
  TCut goodfit("dequal.TrkQualDeM>0.4");

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
  _tn->Project("allrec","sqrt(demcgen.momx^2+demcgen.momy^2+demcgen.momz^2)");
  _tn->Project("tqrec","sqrt(demcgen.momx^2+demcgen.momy^2+demcgen.momz^2)",goodfit);
  _tn->Project("tprtrig","sqrt(demcgen.momx^2+demcgen.momy^2+demcgen.momz^2)",goodfit&&goodtpr);
  _tn->Project("cprtrig","sqrt(demcgen.momx^2+demcgen.momy^2+demcgen.momz^2)",goodfit&&goodcpr);
  _tn->Project("trktrig","sqrt(demcgen.momx^2+demcgen.momy^2+demcgen.momz^2)",goodfit&&goodtrk);
  _tn->Project("alltrig","sqrt(demcgen.momx^2+demcgen.momy^2+demcgen.momz^2)",goodfit&&goodtrig);
  _tn->Project("cctrig","sqrt(demcgen.momx^2+demcgen.momy^2+demcgen.momz^2)",goodfit&&cc);

  // scale by the absolute normalization
  double scalefac =100.0/norm;
  allrec->Scale(scalefac);
  tqrec->Scale(scalefac);
  tprtrig->Scale(scalefac);
  cprtrig->Scale(scalefac);
  trktrig->Scale(scalefac);
  alltrig->Scale(scalefac);
  cctrig->Scale(scalefac);
  _effcan = new TCanvas("effcan","Efficiency",600,600);
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

void TrkAnaPlots::PlotIPA() {
  gStyle->SetOptFit(111111);
  gStyle->SetOptStat("oumr");
  TH1F* trkqual = new TH1F("trkqual","TrkQual",103,-0.01,1.01);
  TH1F* mom = new TH1F("mom","Reco Momentum;P_{reco} (MeV/c)",100,40,56);
  TH1F* momres = new TH1F("momres","Momentum Resolution;P_{reco} - P_{MC} (MeV/c)",100,-5.0,5.0);
  TH1F* nactive = new TH1F("nactive","N Active Straw Hits",121,-0.5,120.5);
  trkqual->SetStats(0);
  _tn->Project("trkqual","dequal.TrkQualDeM");
  _tn->Project("mom","deent.mom","dequal.TrkQualDeM>0.4");
  _tn->Project("momres","deent.mom-sqrt(demcent.momx^2+demcent.momy^2+demcent.momz^2)","dequal.TrkQualDeM>0.4");
  _tn->Project("nactive","de.nactive","dequal.TrkQualDeM>0.4");
  _ipacan = new TCanvas("ipacan","ipacan",800,800);
  _ipacan->Divide(2,2);
  _ipacan->cd(1);
  trkqual->Draw();
  _ipacan->cd(2);
  nactive->Draw();
  _ipacan->cd(3);
  mom->Draw();
  _ipacan->cd(4);
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

void TrkAnaPlots::Upstream() {
  TCut trueup("de.pdg*demc.pdg>0 && demcxit.momz<0");
  TCut uetch("uetch.active");
  TCut trueutch("uetchmc.prel>=0");
  TCut truee("abs(demc.pdg)==11");
  TCut truemu("abs(demc.pdg)==13");
  TH1F* tchmcrel = new TH1F("tchmcrel","TrkCaloHit MC Relation;Calo WRT Track Relationship",8,-1.5,6.5);
  tchmcrel->GetXaxis()->SetBinLabel(1,"none");
  tchmcrel->GetXaxis()->SetBinLabel(2,"same");
  tchmcrel->GetXaxis()->SetBinLabel(3,"daughter");
  tchmcrel->GetXaxis()->SetBinLabel(4,"mother");
  tchmcrel->GetXaxis()->SetBinLabel(5,"sibling");
  tchmcrel->GetXaxis()->SetBinLabel(6,"u-daughter");
  tchmcrel->GetXaxis()->SetBinLabel(7,"u-mother");
  tchmcrel->GetXaxis()->SetBinLabel(8,"u-sibling");
  tchmcrel->SetStats(0);
  TH1F* eutcha = new TH1F("eutcha","Upstream TrkCaloHit",2,-0.5,1.5);
  TH1F* muutcha = new TH1F("muutcha","Upstream TrkCaloHit",2,-0.5,1.5);
  eutcha->GetXaxis()->SetBinLabel(1,"None/Inactive");
  eutcha->GetXaxis()->SetBinLabel(2,"Active");
  eutcha->SetStats(0);
  eutcha->SetLineColor(kBlue);
  muutcha->GetXaxis()->SetBinLabel(1,"None/Inactive");
  muutcha->GetXaxis()->SetBinLabel(2,"Active");
  muutcha->SetStats(0);
  muutcha->SetLineColor(kBlack);
  TH1F* updg = new TH1F("updg","True Upstream Particle PDG code",27,-13.5,13.5);
  updg->GetXaxis()->SetBinLabel(1,"#mu^{+}");
  updg->GetXaxis()->SetBinLabel(3,"e^{+}");
  updg->GetXaxis()->SetBinLabel(27,"#mu^{-}");
  updg->GetXaxis()->SetBinLabel(25,"e^{-}");
  updg->SetStats(0);

  TH1F* eutime = new TH1F("eutime","Upstream fit calo - track time;T_{calo}-T_{track} (ns)",200,-10,10);
  TH1F* muutime = new TH1F("muutime","Upstream fit calo - track time;T_{calo}-T_{track} (ns)",200,-10,10);
  muutime->SetStats(0);
  eutime->SetLineColor(kBlue);
  muutime->SetLineColor(kBlack);
  TH2F* ueevsp = new TH2F("ueevsp","Upstream electron E vs P;Fit mom (MeV/c);CaloCluster EDep (MeV)",25,60,200,25,0,700);
  TH2F* umuevsp = new TH2F("umuevsp","Upstream muon E vs P;Fit mom (MeV/c);CaloCluster EDep (MeV)",25,60,200,25,0,700);
  ueevsp->SetStats(0);
  umuevsp->SetStats(0);

  _tn->Project("eutcha","uetch.active",trueup&&truee);
  _tn->Project("muutcha","uetch.active",trueup&&truemu);
  _tn->Project("tchmcrel","uetchmc.prel",trueup&&uetch);
  _tn->Project("eutime","uetch.ctime+uetch.clen/200.0-uetch.t0",trueup&&uetch&&trueutch&&truee);
  _tn->Project("updg","demc.pdg",trueup);
  _tn->Project("muutime","uetch.ctime+uetch.clen/200.0-uetch.t0",trueup&&uetch&&trueutch&&truemu);
  _tn->Project("ueevsp","uetch.edep:ue.mom",trueup&&uetch&&trueutch&&truee);
  _tn->Project("umuevsp","uetch.edep:ue.mom",trueup&&uetch&&trueutch&&truemu);

  _tp->Project("+eutcha","uetch.active",trueup&&truee);
  _tp->Project("+muutcha","uetch.active",trueup&&truemu);
  _tp->Project("+tchmcrel","uetchmc.prel",trueup&&uetch);
  _tp->Project("+updg","demc.pdg",trueup);
  _tp->Project("+eutime","uetch.ctime+uetch.clen/200.0-uetch.t0",trueup&&uetch&&trueutch&&truee);
  _tp->Project("+muutime","uetch.ctime+uetch.clen/200.0-uetch.t0",trueup&&uetch&&trueutch&&truemu);
  _tp->Project("+ueevsp","uetch.edep:ue.mom",trueup&&uetch&&trueutch&&truee);
  _tp->Project("+umuevsp","uetch.edep:ue.mom",trueup&&uetch&&trueutch&&truemu);

  _ucan = new TCanvas("ucan","Upstream",800,800);
  _ucan->Divide(2,2);
  _ucan->cd(1);
  updg->Draw();
  _ucan->cd(2);
  muutcha->Draw();
  eutcha->Draw("same");
  TLegend* tchleg = new TLegend(0.6,0.7,0.9,0.9);
  tchleg->AddEntry(eutcha,"True electron track ","l");
  tchleg->AddEntry(muutcha,"True muon track","l");
  tchleg->Draw();
  _ucan->cd(3);
  tchmcrel->Draw(); 
  _ucan->cd(4);
  muutime->Fit("gaus");
  eutime->Fit("gaus","","sames");
  TLegend* tleg = new TLegend(0.1,0.7,0.4,0.9);
  tleg->AddEntry(eutime,"True electron track ","l");
  tleg->AddEntry(muutime,"True muon track","l");
  tleg->Draw();

  _uecan = new TCanvas("uecan","Upstream E vs P",800,800);
  _uecan->Divide(1,2);
  _uecan->cd(1);
  gPad->SetLogz();
  ueevsp->Draw("colorz");
  _uecan->cd(2);
  gPad->SetLogz();
  umuevsp->Draw("colorz");
}

void TrkAnaPlots::PBI(unsigned ngen, int charge) {
  TCut truece("demc.gen==2");
    TCut emall = _eminustrig+_eminustq+_livegate+_rpitch+_opa+_upstream+_CRV+_eminuspid+_eminusrmom+truece;
    TCut epall = _eplustrig+_eplustq+_livegate+_rpitch+_opa+_upstream+_CRV+_epluspid+_eplusrmom+truece;
//  TCut emall = _eminustrig+_CRV+truece;
//  TCut epall = _eplustrig+_CRV+truece;
  TH1F* emeff = new TH1F("emeff","Selection Efficiency for #mu^{-}#rightarrowe^{-} vs PBI;Relative PBI",100,0,3.0);
  TH1F* epeff = new TH1F("epeff","Selection Efficiency for #mu^{-}#rightarrowe^{+} vs PBI;Relative PBI",100,0,3.0);
  TH1F* emeffw = new TH1F("emeffw","Signal Efficiency for #mu^{-}#rightarrowe^{-} vs PBI;Relative PBI",100,0,3.0);
  TH1F* epeffw = new TH1F("epeffw","Signal Efficiency for #mu^{-}#rightarrowe^{+} vs PBI;Relative PBI",100,0,3.0);
  emeff->Sumw2();
  epeff->Sumw2();
  emeffw->Sumw2();
  epeffw->Sumw2();
  TF1* ln = new TF1("ln","[0]*TMath::LogNormal(x,[1],[2],[3])",0,4.0);
  // compute mu given sigma, assuming the mean = 1.0
  double sigma = 0.3814; // From IHEP fit
  double mu = -0.5*sigma*sigma;
  double median = exp(mu);
  ln->SetParameters(3.0*ngen/100,sigma,0.0,median);
  TCanvas* effcan = new TCanvas("effcan","effcan",800,600);
  effcan->Divide(1,2);
  if(charge<0){
    _tn->Project("emeff","evtwt.PBIWeight",emall);
    _tn->Project("emeffw","evtwt.PBIWeight",emall*_pbi);
    double neteff = emeffw->Integral()/ngen;
    cout << "tracks passing cuts = " << emeff->GetEntries() << " net signal efficiency = " <<  neteff << endl;
    emeff->Divide(ln);
    emeff->SetMaximum(neteff+0.1);
    emeff->SetMinimum(neteff-0.1);
    effcan->cd(1);
    emeff->Fit("pol1");
    effcan->cd(2);
    emeffw->Scale(1.0/ngen);
    emeffw->Draw();
// compute the convolution efficiency
//    TF1* line = emeff->GetFunction("pol1");
//    line->SetName("line");
//    TF1* prod = new TF1("prod","[&](double *x, double *p){return ln(x)*line(x); }",0,3.0,0);

  } else {
    _tp->Project("epeff","evtwt.PBIWeight",epall);
    _tp->Project("epeffw","evtwt.PBIWeight",epall*_pbi);
    cout << "tracks passing cuts = " << epeff->GetEntries() << " naive efficiency " << epeff->GetEntries()/ngen
    << " net efficiency = " << epeffw->Integral()/ngen << endl;
    epeff->Divide(ln);
    epeff->SetMaximum(0.2);
    epeff->SetMinimum(0.0);
    epeff->Fit("pol1");
  }
}

void TrkAnaPlots::Trigger() {
  TCut down("demcent.momz>0");
  TH1F* nmom = new TH1F("nmom","Downstream Negative Reco Momentum;Momentum (MeV/c)",100,50,220);
  TH1F* ntmom = new TH1F("ntmom","Downstream Negative Reco Momentum;Momentum (MeV/c)",100,50,220);
  TH1F* pmom = new TH1F("pmom","Downstream Positive Reco Momentum;Momentum (MeV/c)",100,50,220);
  TH1F* ptmom = new TH1F("ptmom","Downstream Positive Reco Momentum;Momentum (MeV/c)",100,50,220);
  TH1F* nmommu = new TH1F("nmommu","Downstream Negative Reco Momentum;Momentum (MeV/c)",100,50,220);
  TH1F* ntmommu = new TH1F("ntmommu","Downstream Negative Reco Momentum;Momentum (MeV/c)",100,50,220);
  TH1F* pmommu = new TH1F("pmommu","Downstream Positive Reco Momentum;Momentum (MeV/c)",100,50,220);
  TH1F* ptmommu = new TH1F("ptmommu","Downstream Positive Reco Momentum;Momentum (MeV/c)",100,50,220);
  nmom->SetStats(0);
  ntmom->SetStats(0);
  pmom->SetStats(0);
  ptmom->SetStats(0);
  nmom->SetLineColor(kBlack);
  ntmom->SetLineColor(kRed);
  pmom->SetLineColor(kBlack);
  ptmom->SetLineColor(kRed);
  nmommu->SetStats(0);
  ntmommu->SetStats(0);
  pmommu->SetStats(0);
  ptmommu->SetStats(0);
  nmommu->SetLineColor(kGreen);
  ntmommu->SetLineColor(kCyan);
  pmommu->SetLineColor(kGreen);
  ptmommu->SetLineColor(kCyan);
  TH1F* nd0 = new TH1F("nd0","Downstream Negative Reco d_{0};d_{0} (mm)",100,-500,500);
  TH1F* ntd0 = new TH1F("ntd0","Downstream Negative Reco d_{0};d_{0} (mm)",100,-500,500);
  TH1F* pd0 = new TH1F("pd0","Downstream Positive Reco d_{0};d_{0} (mm)",100,-500,500);
  TH1F* ptd0 = new TH1F("ptd0","Downstream Positive Reco d_{0};d_{0} (mm)",100,-500,500);
  
  TH1F* nd0mu = new TH1F("nd0mu","Downstream Negative Reco d_{0};d_{0} (mm)",100,-500,500);
  TH1F* ntd0mu = new TH1F("ntd0mu","Downstream Negative Reco d_{0};d_{0} (mm)",100,-500,500);
  TH1F* pd0mu = new TH1F("pd0mu","Downstream Positive Reco d_{0};d_{0} (mm)",100,-500,500);
  TH1F* ptd0mu = new TH1F("ptd0mu","Downstream Positive Reco d_{0};d_{0} (mm)",100,-500,500);
  nd0->SetStats(0);
  ntd0->SetStats(0);
  pd0->SetStats(0);
  ptd0->SetStats(0);
  nd0->SetLineColor(kBlack);
  ntd0->SetLineColor(kRed);
  pd0->SetLineColor(kBlack);
  ptd0->SetLineColor(kRed);

  nd0mu->SetStats(0);
  ntd0mu->SetStats(0);
  pd0mu->SetStats(0);
  ptd0mu->SetStats(0);
  nd0mu->SetLineColor(kGreen);
  ntd0mu->SetLineColor(kCyan);
  pd0mu->SetLineColor(kGreen);
  ptd0mu->SetLineColor(kCyan);

  _tn->Project("nmom","deent.mom",_downstream+_eminus);
  _tn->Project("ntmom","deent.mom",_downstream+_eminus+_eminustrig);
  _tp->Project("pmom","deent.mom",_downstream+_eplus);
  _tp->Project("ptmom","deent.mom",_downstream+_eplus+_eplustrig);

  _tn->Project("nmommu","deent.mom",_downstream+_muminus);
  _tn->Project("ntmommu","deent.mom",_downstream+_muminus+_eminustrig);
  _tp->Project("pmommu","deent.mom",_downstream+_muplus);
  _tp->Project("ptmommu","deent.mom",_downstream+_muplus+_eplustrig);

  _tn->Project("nd0","deent.d0",_downstream+_eminus);
  _tn->Project("ntd0","deent.d0",_downstream+_eminus+_eminustrig);
  _tp->Project("pd0","deent.d0",_downstream+_eplus);
  _tp->Project("ptd0","deent.d0",_downstream+_eplus+_eplustrig);

  _tn->Project("nd0mu","deent.d0",_downstream+_muminus);
  _tn->Project("ntd0mu","deent.d0",_downstream+_muminus+_eminustrig);
  _tp->Project("pd0mu","deent.d0",_downstream+_muplus);
  _tp->Project("ptd0mu","deent.d0",_downstream+_muplus+_eplustrig);

  TLegend* tleg = new TLegend(0.6,0.6,0.9,0.9);
  tleg->AddEntry(nmom,"Reconstructed e","L");
  tleg->AddEntry(ntmom,"Triggered e","L");
  tleg->AddEntry(nmommu,"Reconstructed #mu","L");
  tleg->AddEntry(ntmommu,"Triggered #mu","L");

  TCanvas* tcan = new TCanvas("tcan","tcan",800,800);
  tcan->Divide(2,2);
  tcan->cd(1);
  nmom->Draw();
  tleg->Draw();
  ntmom->Draw("same");
  nmommu->Draw("same");
  ntmommu->Draw("same");
  tcan->cd(2);
  pmommu->Draw();
  ptmommu->Draw("same");
  pmom->Draw("same");
  ptmom->Draw("same");
  tcan->cd(3);
  nd0->Draw();
  ntd0->Draw("same");
  nd0mu->Draw("same");
  ntd0mu->Draw("same");
  tcan->cd(4);
  pd0mu->Draw();
  ptd0mu->Draw("same");
  pd0->Draw("same");
  ptd0->Draw("same");
}

void TrkAnaPlots::Alg() {
  TH1F* tprnmom = new TH1F("tprnmom","e^{-} Reco Momentum;Momentum (MeV/c)",100,50,220);
  TH1F* cprnmom = new TH1F("cprnmom","e^{-} Reco Momentum;Momentum (MeV/c)",100,50,220);
  THStack* algnmom = new THStack("algnmom","e^{-} Reco Momentum;Momentum (MeV/c)");
  algnmom->Add(cprnmom);
  algnmom->Add(tprnmom);
  tprnmom->SetFillColor(kRed);
  cprnmom->SetFillColor(kCyan);
  TH1F* tprnmommu = new TH1F("tprnmommu","#mu^{-} Reco Momentum;Momentum (MeV/c)",100,50,220);
  TH1F* cprnmommu = new TH1F("cprnmommu","#mu^{-} Reco Momentum;Momentum (MeV/c)",100,50,220);
  THStack* algnmommu = new THStack("algnmommu","#mu^{-} Reco Momentum;Momentum (MeV/c)");
  algnmommu->Add(cprnmommu);
  algnmommu->Add(tprnmommu);
  tprnmommu->SetFillColor(kRed);
  cprnmommu->SetFillColor(kCyan);
  TH1F* tprpmom = new TH1F("tprpmom","e^{+} Reco Momentum;Momentum (MeV/c)",100,50,220);
  TH1F* cprpmom = new TH1F("cprpmom","e^{+} Reco Momentum;Momentum (MeV/c)",100,50,220);
  THStack* algpmom = new THStack("algpmom","e^{+} Reco Momentum;Momentum (MeV/c)");
  algpmom->Add(cprpmom);
  algpmom->Add(tprpmom);
  tprpmom->SetFillColor(kRed);
  cprpmom->SetFillColor(kCyan);
  TH1F* tprpmommu = new TH1F("tprpmommu","#mu^{+} Reco Momentum;Momentum (MeV/c)",100,50,220);
  TH1F* cprpmommu = new TH1F("cprpmommu","#mu^{+} Reco Momentum;Momentum (MeV/c)",100,50,220);
  THStack* algpmommu = new THStack("algpmommu","#mu^{+} Reco Momentum;Momentum (MeV/c)");
  algpmommu->Add(cprpmommu);
  algpmommu->Add(tprpmommu);
  tprpmommu->SetFillColor(kRed);
  cprpmommu->SetFillColor(kCyan);

  TH1F* tprnd0 = new TH1F("tprnd0","e^{-} Reco d_{0};d_{0} (mm)",100,-500,500);
  TH1F* cprnd0 = new TH1F("cprnd0","e^{-} Reco d_{0};d_{0} (mm)",100,-500,500);
  THStack* algnd0 = new THStack("algnd0","e^{-} Reco d_{0};d_{0} (mm)");
  algnd0->Add(cprnd0);
  algnd0->Add(tprnd0);
  tprnd0->SetFillColor(kRed);
  cprnd0->SetFillColor(kCyan);
  TH1F* tprnd0mu = new TH1F("tprnd0mu","#mu^{-} Reco d_{0};d_{0} (mm)",100,-500,500);
  TH1F* cprnd0mu = new TH1F("cprnd0mu","#mu^{-} Reco d_{0};d_{0} (mm)",100,-500,500);
  THStack* algnd0mu = new THStack("algnd0mu","#mu^{-} Reco d_{0};d_{0} (mm)");
  algnd0mu->Add(cprnd0mu);
  algnd0mu->Add(tprnd0mu);
  tprnd0mu->SetFillColor(kRed);
  cprnd0mu->SetFillColor(kCyan);
  TH1F* tprpd0 = new TH1F("tprpd0","e^{+} Reco d_{0};d_{0} (mm)",100,-500,500);
  TH1F* cprpd0 = new TH1F("cprpd0","e^{+} Reco d_{0};d_{0} (mm)",100,-500,500);
  THStack* algpd0 = new THStack("algpd0","e^{+} Reco d_{0};d_{0} (mm)");
  algpd0->Add(cprpd0);
  algpd0->Add(tprpd0);
  tprpd0->SetFillColor(kRed);
  cprpd0->SetFillColor(kCyan);
  TH1F* tprpd0mu = new TH1F("tprpd0mu","#mu^{+} Reco d_{0};d_{0} (mm)",100,-500,500);
  TH1F* cprpd0mu = new TH1F("cprpd0mu","#mu^{+} Reco d_{0};d_{0} (mm)",100,-500,500);
  THStack* algpd0mu = new THStack("algpd0mu","#mu^{+} Reco d_{0};d_{0} (mm)");
  algpd0mu->Add(cprpd0mu);
  algpd0mu->Add(tprpd0mu);
  tprpd0mu->SetFillColor(kRed);
  cprpd0mu->SetFillColor(kCyan);

  _tn->Project("tprnmom","deent.mom",_eminus+_downstream+_TPR);
  _tn->Project("cprnmom","deent.mom",_eminus+_downstream+_CPR);
  _tp->Project("tprpmom","deent.mom",_eplus+_downstream+_TPR);
  _tp->Project("cprpmom","deent.mom",_eplus+_downstream+_CPR);
  _tn->Project("tprnmommu","deent.mom",_muminus+_downstream+_TPR);
  _tn->Project("cprnmommu","deent.mom",_muminus+_downstream+_CPR);
  _tp->Project("tprpmommu","deent.mom",_muplus+_downstream+_TPR);
  _tp->Project("cprpmommu","deent.mom",_muplus+_downstream+_CPR);

  _tn->Project("tprnd0","deent.d0",_eminus+_downstream+_TPR);
  _tn->Project("cprnd0","deent.d0",_eminus+_downstream+_CPR);
  _tp->Project("tprpd0","deent.d0",_eplus+_downstream+_TPR);
  _tp->Project("cprpd0","deent.d0",_eplus+_downstream+_CPR);
  _tn->Project("tprnd0mu","deent.d0",_muminus+_downstream+_TPR);
  _tn->Project("cprnd0mu","deent.d0",_muminus+_downstream+_CPR);
  _tp->Project("tprpd0mu","deent.d0",_muplus+_downstream+_TPR);
  _tp->Project("cprpd0mu","deent.d0",_muplus+_downstream+_CPR);

  TLegend* aleg = new TLegend(0.6,0.7,0.9,0.9);
  aleg->AddEntry(tprnmom,"TrkPatRec","F");
  aleg->AddEntry(cprnmom,"CalPatRec","F");
  
  TCanvas* acane = new TCanvas("acane","acane",800,800);
  acane->Divide(2,2);
  acane->cd(1);
  algnmom->Draw();
  aleg->Draw();
  acane->cd(2);
  algpmom->Draw();
  acane->cd(3);
  algnd0->Draw();
  acane->cd(4);
  algpd0->Draw();
  TCanvas* acanmu = new TCanvas("acanmu","acanmu",800,800);
  acanmu->Divide(2,2);
  acanmu->cd(1);
  algnmommu->Draw();
  aleg->Draw();
  acanmu->cd(2);
  algpmommu->Draw();
  acanmu->cd(3);
  algnd0mu->Draw();
  acanmu->cd(4);
  algpd0mu->Draw();
}

void TrkAnaPlots::dEdx() {
  TCut momcut("abs(deent.mom-105)<10");
  TH1F* ededx = new TH1F("ededx","Median Active Hit dE/dx;dE/dx (KeV/mm)",100,0,0.4);
  TH1F* mudedx = new TH1F("mudedx","Median Active Hit dE/dx;dE/dx (KeV/mm)",100,0,0.4);
  ededx->SetLineColor(kRed);
  mudedx->SetLineColor(kBlue);
  ededx->SetStats(0);
  mudedx->SetStats(0);
  _tn->Project("ededx","1000*detrkpid.mdedx",_eminus+_downstream+momcut);
  _tn->Project("mudedx","1000*detrkpid.mdedx",_muminus+_downstream+momcut);
  _tp->Project("+ededx","1000*detrkpid.mdedx",_eplus+_downstream+momcut);
  _tp->Project("+mudedx","1000*detrkpid.mdedx",_muplus+_downstream+momcut);
  TCanvas* dedxcan = new TCanvas("dedxcan","dedxcan",600,600);
  ededx->Draw();
  mudedx->Draw("same");
  TLegend* dedxleg = new TLegend(0.6,0.7,0.9,0.9);
  dedxleg->AddEntry(ededx,"True Electrons","L");
  dedxleg->AddEntry(mudedx,"True Muons","L");
  dedxleg->Draw();
}
