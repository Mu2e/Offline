#include "TH1F.h"
#include "TF1.h"
#include "TTree.h"
#include "TLegend.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TCut.h"
#include "TMath.h"
#include "TProfile.h"
#include "TDirectory.h"

void SDTest(TTree* sddiag, const char* page ="adc",unsigned NADC=12,TCut cut=TCut()) {
  TString spage(page);
  if(spage == "adc") {
    TCanvas* adc = new TCanvas("adc","ADC",800,600);
    TH2F* adcwf = new TH2F("adcwf","ADC Waveform;Digitization;Value",NADC,-0.5,NADC-0.5,1024,1300,4095.5);
    TProfile* adcwfp = new TProfile("adcwfp","ADC Waveform;Digitization;Value",NADC,-0.5,NADC-0.5,1300,4095.5,"S");
    adcwf->SetStats(0);
    adcwfp->SetStats(0);
    char name[15];
    for(unsigned iadc=0;iadc<NADC;++iadc){
      snprintf(name,15,"adc[%u]:%u",iadc,iadc);
      std::cout << "name = " << name << std::endl;
      sddiag->Project("+adcwf",name,cut);
      sddiag->Project("+adcwfp",name,cut);
    }
    adc->Divide(1,2);
    adc->cd(1);
    adcwf->Draw("box");
    adc->cd(2);
    adcwfp->Draw();
  } else if(spage == "tdc"){
    TCanvas* tdc = new TCanvas("tdc","TDC",800,600);
    TH1F* dtdc = new TH1F("dtdc","#Delta TDC",200,-500,500);
    sddiag->Project("dtdc","tdccal-tdchv");
    dtdc->Draw();
  } else if(spage=="dvdt") {
    TH2F* dvdt = new TH2F("dvdt","#Delta V vs #Delta t;#Delta t (nsec);#Delta V (mm)",100,-6,6,100,-1200,1200);
    sddiag->Project("dvdt","wdisthv-wdistcal:xtimehv-xtimecal","vcrosscal>0&&vcrosshv>0");
    dvdt->FitSlicesY(0,20,80);
    TProfile* dtdcdt = new TProfile("dtdcdt","#Delta TDC vs #Delta t;#Delta t (ns);#Delta TDC",100,-6,6,-200,200);
    sddiag->Project("dvdt","wdisthv-wdistcal:xtimehv-xtimecal","vcrosscal>0&&vcrosshv>0");
    sddiag->Project("dtdcdt","tdchv-tdccal:xtimehv-xtimecal","vcrosscal>0&&vcrosshv>0");
    dvdt->FitSlicesY(0,20,80);
    TCanvas* dvdtcan = new TCanvas("dvdtcan","dvdtcan",800,600);
    dvdtcan->Divide(2,2);
    dvdtcan->cd(1);
    dvdt->Draw("box");
    TH1D* dvdt_1 = (TH1D*)gDirectory->Get("dvdt_1");
    dvdt_1->SetTitle("Average #Delta V vs #Delta t;#Delta t (ns);#Delta V (mm)");
    dvdtcan->cd(2);
    dvdt_1->Fit("pol1");
    dvdtcan->cd(3);
    dtdcdt->Fit("pol1");

  } else if(spage=="clusters") {

    TH1F* nclust = new TH1F("nclust","N Cluster in Straw",101,-.05,100.5);
    TH1F* iclust = new TH1F("iclust","Threshold - Earliest Cluster Index;I_{thresh}-I_{earliest}",21,-.05,20.5);
    TH1F* tclust = new TH1F("tclust","Threshold - Earliest Cluster Time;T_{thresh}-T_{earliest}(ns)",100,0.0,30.0);
    TH2F* etvsd = new TH2F("etvsd","Earliest Cluster Time vs MC DOCA;MC DOCA (mm);T_{earliest}-T_{MC}",50,0,2.5,50,0,45.0);
    TH2F* ttvsd = new TH2F("ttvsd","Threshold Cluster Time vs MC DOCA;MC DOCA (mm);T_{thresh}-T_{MC}",50,0,2.5,50,0,45.0);
    TH2F* xtvsd = new TH2F("xtvsd","Threhold Xing Drift vs MC DOCA;MC DOCA (mm);T_{xing}-T_{MC}",50,0,2.5,50,0,45.0);

    sddiag->Project("iclust","iclustcal","mcmom>90");
    sddiag->Project("nclust","nclustcal","mcmom>90");
    sddiag->Project("tclust","tctimecal-ectimecal","mcmom>90");
    sddiag->Project("etvsd","ectimecal-mctime%1695:mcdca","mcmom>50&ectime>300");
    sddiag->Project("ttvsd","tctimecal-mctime%1695:mcdca","mcmom>50&ectime>300");
    sddiag->Project("xtvsd","xtimecal-mctime%1695:mcdca","mcmom>50&ectime>300");

    TCanvas* ccan = new TCanvas("ccan","ccan",800,600);
    ccan->Divide(3,2);
    ccan->cd(1);
    nclust->Draw();
    ccan->cd(2);
    iclust->Draw();
    ccan->cd(3);
    tclust->Draw();
    ccan->cd(4);
    etvsd->Draw("box");
    ccan->cd(5);
    ttvsd->Draw("box");
    ccan->cd(6);
    xtvsd->Draw("box");
  } else if(spage == "ionize") {
    TH1F* dke = new TH1F("dke","Compton and #delta-Ray e Kinetic Energy;Energy (KeV)",100,0,10.0);
    TH1F* dbg = new TH1F("dbg","Compton and #delta-Ray e #beta#times#gamma;#beta#times#gamma",100,0,1.0);
    TH1F* pbg = new TH1F("pbg","Proton #beta#times#gamma;#beta#times#gamma",100,0,1.0);

    sddiag->Project("dke","(sqrt(mcmom^2+0.511^2)-0.511)*1000.0","mcpdg==11&&mcproc<20");
    sddiag->Project("dbg","mcmom/0.511","mcpdg==11&&mcproc<20");
    sddiag->Project("pbg","mcmom/938.27","mcpdg==2212");
    TCanvas* ican = new TCanvas("ican","ican",800,600);
    ican->Divide(2,2);
    ican->cd(1);
    dke->Draw();
    ican->cd(2);
    dbg->Draw();
    ican->cd(3);
    pbg->Draw();
  }
}
