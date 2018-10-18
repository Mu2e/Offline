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
#include <string>
using std::string;
class TimeClusterDiag  {
  public:
    TimeClusterDiag(TTree* tcdiag) : _tcdiag(tcdiag),
    _goodCENHits("ceclust.nce > 14"),
    _goodCETime("ceclust.time > 500.0"),
    _goodReco("besttc.nhits >= 10 "),
    _goodCEReco("besttc.ncehits > 0 "),
    _goodCalo("besttc.ecalo>50.0"),
    _disk1("besttc.cogz<2000.0"),
    _disk2("besttc.cogz>2000.0"),
    _cehit("tchinfo._mcrel==0&&tchinfo._mcgen==2"),
    _effcan(0), _btccan(0), _cecan(0), _tcan(0), _ctcan(0)
  {
    _goodCE = _goodCENHits+_goodCETime;
  }

    TTree* _tcdiag;
    TCut _goodCENHits, _goodCETime, _goodCE;
    TCut _goodReco, _goodCEReco;
    TCut _goodCalo, _disk1, _disk2;
    TCut _cehit;
    TCanvas *_effcan, *_btccan, *_cecan, *_tcan, *_ctcan, *_hcan;

    void Efficiency();
    void BestTC();
    void CE();
    void time();
    void ctime();
    void HitTimeRes();
    void save(const char* suffix=".png");
};

void TimeClusterDiag::Efficiency() {
  TH1F* tct = new TH1F("tct","Time Cluster Reco Efficiency vs Time;Time (ns); Efficiency",100,0.0,1695.0);
  TH1F* tctr = new TH1F("tctr","Time Cluster Reco Efficiency vs Time;Time (ns); Efficiency",100,0.0,1695.0);
  TH1F* tcte = new TH1F("tcte","Time Cluster Reco Efficiency vs Time;Time (ns); Efficiency",100,0.0,1695.0);
  TH1F* tcn = new TH1F("tcn","NHits Cluster Reco Efficiency vs NHits;NHits; Efficiency",100,-0.5,99.5);
  TH1F* tcnr = new TH1F("tcnr","NHits Cluster Reco Efficiency vs NHits;NHits; Efficiency",100,-0.5,99.5);
  TH1F* tcne = new TH1F("tcne","NHits Cluster Reco Efficiency vs NHits;NHits; Efficiency",100,-0.5,99.5);

  tcte->Sumw2();
  tcne->Sumw2();
  _tcdiag->Project("tct","ceclust.time",_goodCENHits);
  _tcdiag->Project("tctr","ceclust.time",_goodCENHits+_goodCEReco);
  for(int ibin =0;ibin < tctr->GetNbinsX();++ibin){
    if(tct->GetBinContent(ibin) > 0){
      double nr = tctr->GetBinContent(ibin);
      double nu = tct->GetBinContent(ibin) - nr;
      tcte->SetBinContent(ibin,nr/(nu+nr));
      tcte->SetBinError(ibin,sqrt(nu*nr/pow(nu+nr,3)));
    }
  }
  _tcdiag->Project("tcn","ceclust.nce",_goodCETime);
  _tcdiag->Project("tcnr","ceclust.nce",_goodCETime+_goodCEReco);
  for(int ibin =0;ibin < tcnr->GetNbinsX();++ibin){
    if(tcn->GetBinContent(ibin) > 0){
      double nr = tcnr->GetBinContent(ibin);
      double nu = tcn->GetBinContent(ibin) - nr;
      tcne->SetBinContent(ibin,nr/(nu+nr));
      tcne->SetBinError(ibin,sqrt(nu*nr/pow(nu+nr,3)));
    }
  }
  _effcan = new TCanvas("effcan","Efficiency",800,600);
  _effcan->Divide(2,1);
  _effcan->cd(1);
  tcte->Draw();
  _effcan->cd(2);
  tcne->Draw();


}

void TimeClusterDiag::BestTC(){
  TH1F* nhit = new TH1F("nhit","BestTC N hits",100,-0.5,99.5);
  TH1F* nhitc = new TH1F("nhitc","BestTC N hits",100,-0.5,99.5);
//  TH1F* ncehit = new TH1F("ncehit","BestTC N CE hits",100,-0.5,99.5);
  TH1F* hpur = new TH1F("hpur","BestTC CE Hit Purity",100,0,1.0001);
  TH1F* hpurc = new TH1F("hpurc","BestTC CE Hit Purity",100,0,1.0001);
//
  TH1F* heff = new TH1F("heff","BestTC CE Hit Efficiency",100,0,1.0001);
  TH1F* heffc = new TH1F("heffc","BestTC CE Hit Efficiency",100,0,1.0001);
//  TH1F* ecalo = new TH1F("ecalo","BestTC Calo Cluster Energy;E_{Calo} (MeV)",100,0,120.0);
//  TH1F* time = new TH1F("time","BestTC Time;Average Cluster Time (ns)",100,0.0,1695.0);
  TH1F* tres = new TH1F("tres","BestTC Time Resolution;Cluster Time - Ce Time (ns)",100,-20,20.0);
  TH1F* tresc = new TH1F("tresc","BestTC Time Resolution;Cluster Time - Ce Time (ns)",100,-20,20.0);
  nhit->SetLineColor(kBlue);
  nhitc->SetLineColor(kGreen);
  hpur->SetLineColor(kBlue);
  hpurc->SetLineColor(kGreen);
  heff->SetLineColor(kBlue);
  heffc->SetLineColor(kGreen);
  tres->SetLineColor(kBlue);
  tresc->SetLineColor(kGreen);

//  _tcdiag->Project("time","besttc.time",_goodReco&&_goodCENHits&&_goodCEReco);
  _tcdiag->Project("tres","besttc.time-fmod(mcmidt0,1695)",_goodReco&&_goodCENHits&&_goodCEReco);
  _tcdiag->Project("tresc","besttc.time-fmod(mcmidt0,1695)",_goodReco&&_goodCENHits&&_goodCEReco&&_goodCalo);
  _tcdiag->Project("nhit","besttc.nhits",_goodReco&&_goodCENHits&&_goodCEReco);
  _tcdiag->Project("nhitc","besttc.nhits",_goodReco&&_goodCENHits&&_goodCEReco&&_goodCalo);
//  _tcdiag->Project("ncehit","besttc.ncehits",_goodReco&&_goodCENHits&&_goodCEReco);
  _tcdiag->Project("hpur","besttc.ncehits/besttc.nhits",_goodReco&&_goodCENHits&&_goodCEReco);
  _tcdiag->Project("hpurc","besttc.ncehits/besttc.nhits",_goodReco&&_goodCENHits&&_goodCEReco&&_goodCalo);
//
  _tcdiag->Project("heff","besttc.ncehits/ceclust.nhits",_goodReco&&_goodCENHits&&_goodCEReco);
  _tcdiag->Project("heffc","besttc.ncehits/ceclust.nhits",_goodReco&&_goodCENHits&&_goodCEReco&&_goodCalo);
//  _tcdiag->Project("ecalo","besttc.ecalo",_goodReco&&_goodCENHits&&_goodCEReco);
  _btccan = new TCanvas("btccan","BestTC",800,800);
  TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
  leg->AddEntry(nhit,"All TC","L");
  leg->AddEntry(nhitc,"Calo TC","L");
  _btccan->Divide(2,2);
  _btccan->cd(1);
  nhit->Draw();
  nhitc->Draw("same");
  leg->Draw();
  _btccan->cd(2);
  hpur->Draw();
  hpurc->Draw("same");
  _btccan->cd(3);
  heff->Draw();
  heffc->Draw("same");
  _btccan->cd(4);
  tres->Draw();
  tresc->Draw("same");

}

void TimeClusterDiag::CE(){
  TH1F* nhitce = new TH1F("nhitce","N CE hits",100,-0.5,99.5);
  TH1F* timece = new TH1F("timece","CE Time;Average hit time (ns)",100,0.0,1695.0);
  TH1F* dphice = new TH1F("dphice","#Phi extent of CE hits;#phi extent (rad)",100,0.0,6.3);
  TH1F* rmaxce = new TH1F("rmaxce","Max #Rho of CE hits;#rho max (mm)",100,300.0,800.0);
  TH1F* rmince = new TH1F("rmince","Min #Rho of CE hits;#rho max (mm)",100,300.0,800.0);
  _tcdiag->Project("timece","ceclust.time");
  _tcdiag->Project("nhitce","ceclust.nce");
  _tcdiag->Project("dphice","ceclust.maxdphi",_goodCE);
  _tcdiag->Project("rmaxce","ceclust.maxrho",_goodCE);
  _tcdiag->Project("rmince","ceclust.minrho",_goodCE);
  _cecan = new TCanvas("cecan","CE",1000,800);
  _cecan->Divide(3,2);
  _cecan->cd(1);
  nhitce->Draw();
  _cecan->cd(2);
  timece->Draw();
  _cecan->cd(3);
  dphice->Draw();
  _cecan->cd(4);
  rmince->Draw();
  _cecan->cd(5);
  rmaxce->Draw();
}

void TimeClusterDiag::time() {
  TH1F* calodt1 = new TH1F("calodt1","Calo Cluster time - CE time, Disk 1;#Delta t (ns)",100,-15.0,15.0);
  TH1F* calodt2 = new TH1F("calodt2","Calo Cluster time - CE time, Disk 2;#Delta t (ns)",100,-15.0,15.0);
  _tcdiag->Project("calodt1","besttc.tcalo-mcmidt0",_goodCE+_goodReco+_goodCalo+_disk1);
  _tcdiag->Project("calodt2","besttc.tcalo-mcmidt0",_goodCE+_goodReco+_goodCalo+_disk2);

  calodt1->Fit("gaus","qO");
  calodt2->Fit("gaus","qO");
  
  TProfile* hdtzp = new TProfile("hdtzp","Hit time - CE time vs z;z (mm);#Delta t (ns)",50,-1600,1600,-25,75);
  hdtzp->SetStats(1);
  TH2F* hdtz = new TH2F("hdtz","Hit time - CE time vs z;z (mm);#Delta t (ns)",50,-1600,1600,50,-25,75);
  hdtz->SetStats(0);
  _tcdiag->Project("hdtz","tchinfo._time - mcmidt0:tchinfo._z",_goodCE+_goodReco+_cehit);
  _tcdiag->Project("hdtzp","tchinfo._time - mcmidt0:tchinfo._z",_goodCE+_goodReco+_cehit);
  hdtzp->SetMarkerStyle(20);
  hdtzp->SetMarkerColor(kBlack);

  TH1F* htres = new TH1F("htres","Hit time - CE time;#Delta t (ns)",100,-30,30);
  TH1F* chtres = new TH1F("chtres","Hit time - CE time;#Delta t (ns)",100,-30,30);
  htres->SetLineColor(kGreen);
//  htres->SetStats(0);
  chtres->SetLineColor(kBlue);
  _tcdiag->Project("htres","tchinfo._time -22.5 - mcmidt0",_goodCE+_goodReco+_cehit);
  _tcdiag->Project("chtres","tchinfo._time -22.5 - tchinfo._z*0.0047 -mcmidt0",_goodCE+_goodReco+_cehit);

  _tcan = new TCanvas("tcan","Time",800,800);
  _tcan->Divide(2,2);
  _tcan->cd(1);
  calodt1->Draw();
  _tcan->cd(2);
  calodt2->Draw();
  _tcan->cd(3);
  hdtz->Draw("colorz");
  hdtzp->Fit("pol1","+","sames");
  _tcan->cd(4);
  chtres->Draw();
  htres->Draw("sames");
  TLegend* hleg = new TLegend(0.5,0.2,0.9,0.5);
  hleg->AddEntry(htres,"Straw Hit Times","L");
  hleg->AddEntry(chtres,"Corrected Straw Hit Times","L");
  hleg->Draw();

}

void TimeClusterDiag::ctime() {
  TH1F* tcdt = new TH1F("tcdt","Time Cluster time - CE time;#Delta t (ns)",100,-20,20);
  TH1F* tcdtc = new TH1F("tcdtc","Time Cluster time - CE time;#Delta t (ns)",100,-20,20);
  tcdt->SetLineColor(kGreen);
  tcdtc->SetLineColor(kBlue);
  _tcdiag->Project("tcdt","besttc.time-ceclust.time",_goodCE+_goodReco+(!_goodCalo));
  _tcdiag->Project("tcdtc","besttc.time-ceclust.time",_goodCE+_goodReco+_goodCalo);

  TH1F* tcdtp = new TH1F("tcdtp","Time Cluster time - CE time pull",100,-20,20);
  TH1F* tcdtpc = new TH1F("tcdtpc","Time Cluster time - CE timepull",100,-20,20);
  tcdtp->SetLineColor(kGreen);
  tcdtpc->SetLineColor(kBlue);
  _tcdiag->Project("tcdtp","(besttc.time-ceclust.time)/besttc.terr",_goodCE+_goodReco+(!_goodCalo));
  _tcdiag->Project("tcdtpc","(besttc.time-ceclust.time)/besttc.terr",_goodCE+_goodReco+_goodCalo);

  _ctcan = new TCanvas("ctcan","Cluster Time",800,800);
  _ctcan->Divide(2,2);
  _ctcan->cd(1);
  tcdtc->Fit("gaus");
  _ctcan->cd(2);
  tcdt->Fit("gaus");
  TLegend* tcleg = new TLegend(0.1,0.5,0.4,0.9);
  tcleg->AddEntry(tcdtc,"Trk Hits + Calo Cluster","L");
  tcleg->AddEntry(tcdt,"Trk Hits only","L");
  tcleg->Draw();
  _ctcan->cd(3);
  tcdtpc->Fit("gaus");
  _ctcan->cd(4);
  tcdtp->Fit("gaus");
}

void TimeClusterDiag::HitTimeRes() {
  TH1F* raw = new TH1F("raw","Combo Hit T0 Resolution;T0 (nsec)",100,-30.0,30);
  TH1F* flt = new TH1F("flt","Combo Hit T0 Resolution;T0 (nsec)",100,-30.0,30);
  TH1F* TOT = new TH1F("TOT","Combo Hit T0 Resolution;T0 (nsec)",100,-30.0,30);
  raw->SetLineColor(kRed);
  flt->SetLineColor(kGreen);
  TOT->SetLineColor(kBlue);
  raw->SetStats(0);
  flt->SetStats(0);
  TOT->SetStats(0);
  _tcdiag->Project("raw","tchinfo._time - mcmidt0 - 22.5",_goodCE+_goodReco+_goodCEReco+_cehit);
  _tcdiag->Project("flt","tchinfo._time -22.6 - tchinfo._z*0.00535 -mcmidt0",_goodCE+_goodReco+_goodCEReco+_cehit);
  _tcdiag->Project("TOT","tchinfo._dt+besttc.time - mcmidt0",_goodCE+_goodReco+_goodCEReco+_cehit);
  _hcan = new TCanvas("hcan","hcan",600,600);
  _hcan->Divide(1,1);
  _hcan->cd(1);
  TOT->Draw();
  raw->Draw("same");
  flt->Draw("same");
  TLegend* leg = new TLegend(0.6,0.7,0.9,0.9);
  char cap[40];
  snprintf(cap,40,"Raw Hit Time, RMS=%3.2f",raw->GetRMS());
  leg->AddEntry(raw,cap,"L");
  snprintf(cap,40,"Flt Corr Time, RMS=%3.2f",flt->GetRMS());
  leg->AddEntry(flt,cap,"L");
  snprintf(cap,40,"TOT+Flt Corr Time, RMS=%3.2f",TOT->GetRMS());
  leg->AddEntry(TOT,cap,"L");
  leg->Draw();
}
  
void TimeClusterDiag::save(const char* suffix) {
  string ss(suffix);
  string cfname;
  if(_effcan){
    cfname = string(_effcan->GetName())+ss;
    _effcan->SaveAs(cfname.c_str());
  }
  if(_btccan){
    cfname = string(_btccan->GetName())+ss;
    _btccan->SaveAs(cfname.c_str());
  }
  if(_cecan){
    cfname = string(_cecan->GetName())+ss;
    _cecan->SaveAs(cfname.c_str());
  }
  if(_tcan){
    cfname = string(_tcan->GetName())+ss;
    _tcan->SaveAs(cfname.c_str());
  }
  if(_ctcan){
    cfname = string(_ctcan->GetName())+ss;
    _ctcan->SaveAs(cfname.c_str());
  }
  if(_hcan){
    cfname = string(_hcan->GetName())+ss;
    _hcan->SaveAs(cfname.c_str());
  }
}
