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

class TimeClusterDiag  {
  public:
    TimeClusterDiag(TTree* tcdiag) : _tcdiag(tcdiag),
    _goodCENHits("ceclust.nce > 14"),
    _goodCETime("ceclust.time > 500.0"),
    _goodReco("besttc.nhits > 0 "),
    _goodCalo("besttc.ecalo>50.0"),
    _disk1("besttc.cogz<2000.0"),
    _disk2("besttc.cogz>2000.0"),
    _cehit("tchinfo._mcproc==56")
  {
    _goodCE = _goodCENHits+_goodCETime;
    _goodCEReco = _goodReco +TCut("besttc.ncehits/besttc.nhits>0.5");
  }

    TTree* _tcdiag;
    TCut _goodCENHits, _goodCETime, _goodCE;
    TCut _goodReco, _goodCEReco;
    TCut _goodCalo, _disk1, _disk2;
    TCut _cehit;

    void Efficiency();
    void BestTC();
    void CE();
    void time();
    void ctime();
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
  TCanvas* effcan = new TCanvas("effcan","Efficiency",800,600);
  effcan->Divide(2,1);
  effcan->cd(1);
  tcte->Draw();
  effcan->cd(2);
  tcne->Draw();


}

void TimeClusterDiag::BestTC(){
  TH1F* nhit = new TH1F("nhit","BestTC N hits",100,-0.5,99.5);
  TH1F* ncehit = new TH1F("ncehit","BestTC N CE hits",100,-0.5,99.5);
  TH1F* hpur = new TH1F("hpur","BestTC CE Hit Fraction",100,0,1.0001);

  TH1F* time = new TH1F("time","BestTC Time;Average Cluster Time (ns)",100,0.0,1695.0);
  _tcdiag->Project("time","besttc.time",_goodReco);
  _tcdiag->Project("nhit","besttc.nhits",_goodReco);
  _tcdiag->Project("ncehit","besttc.ncehits",_goodReco);
  _tcdiag->Project("hpur","besttc.ncehits/besttc.nhits",_goodReco);
  TCanvas* btccan = new TCanvas("btccan","BestTC",800,800);
  btccan->Divide(2,2);
  btccan->cd(1);
  nhit->Draw();
  btccan->cd(2);
  hpur->Draw();
  btccan->cd(3);
  time->Draw();
}

void TimeClusterDiag::CE(){
  TH1F* nhitce = new TH1F("nhitce","N CE hits",100,-0.5,99.5);
  TH1F* timece = new TH1F("timece","CE Time;Average hit time (ns)",100,0.0,1695.0);
  TH1F* dphice = new TH1F("dphice","#Phi extent of CE hits;#phi extent (rad)",100,0.0,1.5);
  TH1F* rmaxce = new TH1F("rmaxce","Max #Rho of CE hits;#rho max (mm)",100,300.0,800.0);
  TH1F* rmince = new TH1F("rmince","Min #Rho of CE hits;#rho max (mm)",100,300.0,800.0);
  _tcdiag->Project("timece","ceclust.time");
  _tcdiag->Project("nhitce","ceclust.nce");
  _tcdiag->Project("dphice","ceclust.maxdphi",_goodCE);
  _tcdiag->Project("rmaxce","ceclust.maxrho",_goodCE);
  _tcdiag->Project("rmince","ceclust.minrho",_goodCE);
  TCanvas* cecan = new TCanvas("cecan","CE",1000,800);
  cecan->Divide(3,2);
  cecan->cd(1);
  nhitce->Draw();
  cecan->cd(2);
  timece->Draw();
  cecan->cd(3);
  dphice->Draw();
  cecan->cd(4);
  rmince->Draw();
  cecan->cd(5);
  rmaxce->Draw();
}

void TimeClusterDiag::time() {
  TH1F* calodt1 = new TH1F("calodt1","Calo Cluster time - CE time, Disk 1;#Delta t (ns)",100,0.0,25.0);
  TH1F* calodt2 = new TH1F("calodt2","Calo Cluster time - CE time, Disk 2;#Delta t (ns)",100,0.0,25.0);
  _tcdiag->Project("calodt1","besttc.tcalo-ceclust.time",_goodCE+_goodReco+_goodCalo+_disk1);
  _tcdiag->Project("calodt2","besttc.tcalo-ceclust.time",_goodCE+_goodReco+_goodCalo+_disk2);

  calodt1->Fit("gaus","qO");
  calodt2->Fit("gaus","qO");
  
  TProfile* hdtzp = new TProfile("hdtzp","Hit time - CE time vs z;z (mm);#Delta t (ns)",50,-1600,1600,-25,75);
  hdtzp->SetStats(1);
  TH2F* hdtz = new TH2F("hdtz","Hit time - CE time vs z;z (mm);#Delta t (ns)",50,-1600,1600,50,-25,75);
  hdtz->SetStats(0);
  _tcdiag->Project("hdtz","tchinfo._time - ceclust.time:tchinfo._z",_goodCE+_goodReco+_cehit);
  _tcdiag->Project("hdtzp","tchinfo._time - ceclust.time:tchinfo._z",_goodCE+_goodReco+_cehit);
  hdtzp->SetMarkerStyle(20);
  hdtzp->SetMarkerColor(kBlack);

  TH1F* htres = new TH1F("htres","Hit time - CE time;#Delta t (ns)",100,-30,30);
  TH1F* chtres = new TH1F("chtres","Hit time - CE time;#Delta t (ns)",100,-30,30);
  htres->SetLineColor(kGreen);
//  htres->SetStats(0);
  chtres->SetLineColor(kBlue);
  _tcdiag->Project("htres","tchinfo._time -25.5 -ceclust.time",_goodCE+_goodReco+_cehit);
  _tcdiag->Project("chtres","tchinfo._time -25.5 - tchinfo._z*0.0047 -ceclust.time",_goodCE+_goodReco+_cehit);

  TCanvas* tcan = new TCanvas("tcan","Time",800,800);
  tcan->Divide(2,2);
  tcan->cd(1);
  calodt1->Draw();
  tcan->cd(2);
  calodt2->Draw();
  tcan->cd(3);
  hdtz->Draw("colorz");
  hdtzp->Fit("pol1","+","sames");
  tcan->cd(4);
  chtres->Draw();
  htres->Draw("sames");
  TLegend* hleg = new TLegend(0.5,0.2,0.9,0.5);
  hleg->AddEntry(htres,"Straw Hit Times","L");
  hleg->AddEntry(chtres,"Corrected Straw Hit Times","L");
  hleg->Draw();
}

void TimeClusterDiag::ctime() {
  TH1F* tcdt = new TH1F("tcdt","Time Cluster time - CE time",100,-20,20);
  TH1F* tcdtc = new TH1F("tcdtc","Time Cluster time - CE time",100,-20,20);
  tcdt->SetLineColor(kGreen);
  tcdtc->SetLineColor(kBlue);
  _tcdiag->Project("tcdt","besttc.time-ceclust.time",_goodCE+_goodReco+(!_goodCalo));
  _tcdiag->Project("tcdtc","besttc.time-ceclust.time",_goodCE+_goodReco+_goodCalo);


  TCanvas* ctcan = new TCanvas("ctcan","Cluster Time",800,800);
  tcdtc->Draw();
  tcdt->Draw("same");
  TLegend* tcleg = new TLegend(0.1,0.5,0.4,0.9);
  tcleg->AddEntry(tcdtc,"Trk Hits + Calo Cluster","L");
  tcleg->AddEntry(tcdt,"Trk Hits only","L");
  tcleg->Draw();
}
