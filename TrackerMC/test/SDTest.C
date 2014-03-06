#include "TH1F.h"
#include "TF1.h"
#include "TTree.h"
#include "TLegend.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TCut.h"
#include "TMath.h"
#include "TProfile.h"
#include "TDirectory.h"

void SDTest(TTree* sddiag, char* page ="adc",TCut cut=TCut()) {
  TString spage(page);
  if(spage == "adc") {
    TCanvas* adc = new TCanvas("adc","ADC",800,800);
    TH2F* adcwf = new TH2F("adcwf","ADC Waveform;Digitization;Value",7,-0.5,6.5,1024,-0.5,4095.5);
    TProfile* adcwfp = new TProfile("adcwfp","ADC Waveform;Digitization;Value",7,-0.5,6.5,-0.5,4095.5,"S");
    adcwf->SetStats(0);
    adcwfp->SetStats(0);
    char name[10];
    for(size_t iadc=0;iadc<7;++iadc){
      snprintf(name,10,"adc[%u]:%u",iadc,iadc);
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
    TCanvas* tdc = new TCanvas("tdc","TDC",800,800);
    TH1F* dtdc = new TH1F("dtdc","#Delta TDC",200,-500,500);
    sddiag->Project("dtdc","tdc1-tdc0");
    dtdc->Draw();
  }
}

