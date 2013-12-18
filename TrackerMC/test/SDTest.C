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

void SDTest(TTree* sddiag) {
  TCanvas* adc = new TCanvas("adc","ADC",800,800);
  TH2F* adcwf = new TH2F("adcwf","ADC Waveform;Digitization;Value",10,-0.5,9.5,1024,-0.5,1023.5);
  TProfile* adcwfp = new TProfile("adcwfp","ADC Waveform;Digitization;Value",10,-0.5,9.5,-0.5,1023.5,"S");
  char name[10];
  for(size_t iadc=0;iadc<10;++iadc){
    snprintf(name,10,"adc[%u]:%u",iadc,iadc);
    std::cout << "name = " << name << std::endl;
    sddiag->Project("+adcwf",name);
    sddiag->Project("+adcwfp",name);
  }
  adc->Divide(1,2);
  adc->cd(1);
  adcwf->Draw("box");
  adc->cd(2);
  adcwfp->Draw();

}

