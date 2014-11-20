#include "WLSRunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"

#include "Randomize.hh"

#include "WLSDetectorConstruction.hh"
#include "WLSSteppingAction.hh"
#include "WLSEventAction.hh"

#include <algorithm>
#include <ctime>

#include <TStyle.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TText.h>
#include <TH1D.h>
#include <TH3D.h>
#include <TPaveStats.h>

WLSRunAction::WLSRunAction(int mode) : saveRndm(0), autoSeed(false), _mode(mode)
{
}

WLSRunAction::~WLSRunAction()
{
}

void WLSRunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  G4RunManager::GetRunManager()->SetRandomNumberStoreDir("random/");

  if (autoSeed) {
     // automatic (time-based) random seeds for each run
     G4cout << "*******************" << G4endl;
     G4cout << "*** AUTOSEED ON ***" << G4endl;
     G4cout << "*******************" << G4endl;
     long seeds[2];
     time_t systime = time(NULL);
     seeds[0] = (long) systime;
     seeds[1] = (long) (systime*G4UniformRand());
     CLHEP::HepRandom::setTheSeeds(seeds);
     CLHEP::HepRandom::showEngineStatus();
  } else {
     CLHEP::HepRandom::showEngineStatus();
  }

  if (saveRndm > 0) CLHEP::HepRandom::saveEngineStatus("BeginOfRun.rndm");

}

void WLSRunAction::EndOfRunAction(const G4Run* )
{
  if (saveRndm == 1)
  {
     CLHEP::HepRandom::showEngineStatus();
     CLHEP::HepRandom::saveEngineStatus("endOfRun.rndm");
  }

  if(_mode==0)
  {
    gStyle->SetOptStat(1111);
    TCanvas c1("PEs","PEs",1000,1000);
    c1.Divide(2,2);
    int minbin0[4]={-1, -1, -1, -1};
    int minbin1[4]={-1, -1, -1, -1};
    for(int SiPM=0; SiPM<4; SiPM++)
    {
      int nbins=WLSEventAction::Instance()->GetHistPE(0,SiPM)->GetNbinsX();
      for(int ibin=1; ibin<=nbins; ibin++)
      {
        if(WLSEventAction::Instance()->GetHistPE(0,SiPM)->GetBinContent(ibin)!=0 && minbin0[SiPM]==-1) minbin0[SiPM]=ibin;
        if(WLSEventAction::Instance()->GetHistPE(1,SiPM)->GetBinContent(ibin)!=0 && minbin1[SiPM]==-1) minbin1[SiPM]=ibin;
      }
    }
    int maxbin0[4]={-1, -1, -1, -1};
    int maxbin1[4]={-1, -1, -1, -1};
    for(int SiPM=0; SiPM<4; SiPM++)
    {
      int nbins=WLSEventAction::Instance()->GetHistPE(0,SiPM)->GetNbinsX();
      for(int ibin=nbins; ibin>=1; ibin--)
      {
        if(WLSEventAction::Instance()->GetHistPE(0,SiPM)->GetBinContent(ibin)!=0 && maxbin0[SiPM]==-1) maxbin0[SiPM]=ibin;
        if(WLSEventAction::Instance()->GetHistPE(1,SiPM)->GetBinContent(ibin)!=0 && maxbin1[SiPM]==-1) maxbin1[SiPM]=ibin;
      }
    }
    int minbin=std::min(*std::min_element(minbin0,minbin0+4),*std::min_element(minbin1,minbin1+4));
    int maxbin=std::max(*std::max_element(maxbin0,maxbin0+4),*std::max_element(maxbin1,maxbin1+4));
    for(int SiPM=0; SiPM<4; SiPM++)
    {
      c1.cd(SiPM+1);
      gPad->SetLogy();
      WLSEventAction::Instance()->GetHistPE(0,SiPM)->SetAxisRange(minbin-5,maxbin+5);
      WLSEventAction::Instance()->GetHistPE(0,SiPM)->GetXaxis()->SetTitle("PEs");
      WLSEventAction::Instance()->GetHistPE(0,SiPM)->SetLineColor(1);
      WLSEventAction::Instance()->GetHistPE(0,SiPM)->Draw();
      gPad->Update();
      TPaveStats *stats0 = (TPaveStats*)WLSEventAction::Instance()->GetHistPE(0,SiPM)->FindObject("stats");
      stats0->SetTextColor(1);
      stats0->SetLineColor(1);
      double X1 = stats0->GetX1NDC();
      double Y1 = stats0->GetY1NDC();
      double X2 = stats0->GetX2NDC();
      double Y2 = stats0->GetY2NDC();
      WLSEventAction::Instance()->GetHistPE(1,SiPM)->SetAxisRange(minbin-5,maxbin+5);
      WLSEventAction::Instance()->GetHistPE(1,SiPM)->GetXaxis()->SetTitle("PEs");
      WLSEventAction::Instance()->GetHistPE(1,SiPM)->SetLineColor(2);
      WLSEventAction::Instance()->GetHistPE(1,SiPM)->Draw();
      gPad->Update();
      TPaveStats *stats1 = (TPaveStats*)WLSEventAction::Instance()->GetHistPE(1,SiPM)->FindObject("stats");
      stats1->SetTextColor(2);
      stats1->SetLineColor(2);
      stats1->SetX1NDC(X1);
      stats1->SetY1NDC(Y1-(Y2-Y1));
      stats1->SetX2NDC(X2);
      stats1->SetY2NDC(Y1);
      WLSEventAction::Instance()->GetHistPE(0,SiPM)->Draw();
      WLSEventAction::Instance()->GetHistPE(1,SiPM)->Draw("same");
    }      
    c1.SaveAs("PEs.C");

    TCanvas c2("ArrivalTimes","ArrivalTimes",1000,1000);
    c2.Divide(2,2);
    for(int i=0; i<4; i++)
    {
      minbin0[i]=-1;
      minbin1[i]=-1;
      maxbin0[i]=-1;
      maxbin1[i]=-1;
    }
    for(int SiPM=0; SiPM<4; SiPM++)
    {
      int nbins=WLSEventAction::Instance()->GetHistT(0,SiPM)->GetNbinsX();
      for(int ibin=1; ibin<=nbins; ibin++)
      {
        if(WLSEventAction::Instance()->GetHistT(0,SiPM)->GetBinContent(ibin)!=0 && minbin0[SiPM]==-1) minbin0[SiPM]=ibin;
        if(WLSEventAction::Instance()->GetHistT(1,SiPM)->GetBinContent(ibin)!=0 && minbin1[SiPM]==-1) minbin1[SiPM]=ibin;
      }
    }
    for(int SiPM=0; SiPM<4; SiPM++)
    {
      int nbins=WLSEventAction::Instance()->GetHistT(0,SiPM)->GetNbinsX();
      for(int ibin=nbins; ibin>=1; ibin--)
      {
        if(WLSEventAction::Instance()->GetHistT(0,SiPM)->GetBinContent(ibin)!=0 && maxbin0[SiPM]==-1) maxbin0[SiPM]=ibin;
        if(WLSEventAction::Instance()->GetHistT(1,SiPM)->GetBinContent(ibin)!=0 && maxbin1[SiPM]==-1) maxbin1[SiPM]=ibin;
      }
    }
    minbin=std::min(*std::min_element(minbin0,minbin0+4),*std::min_element(minbin1,minbin1+4));
    maxbin=std::max(*std::max_element(maxbin0,maxbin0+4),*std::max_element(maxbin1,maxbin1+4));
    for(int SiPM=0; SiPM<4; SiPM++)
    {
      c2.cd(SiPM+1);
      gPad->SetLogy();
      WLSEventAction::Instance()->GetHistT(0,SiPM)->SetAxisRange(minbin-5,maxbin+5);
      WLSEventAction::Instance()->GetHistT(0,SiPM)->GetXaxis()->SetTitle("t [ns]");
      WLSEventAction::Instance()->GetHistT(0,SiPM)->SetLineColor(1);
      WLSEventAction::Instance()->GetHistT(0,SiPM)->Draw();
      gPad->Update();
      TPaveStats *stats0 = (TPaveStats*)WLSEventAction::Instance()->GetHistT(0,SiPM)->FindObject("stats");
      stats0->SetTextColor(1);
      stats0->SetLineColor(1);
      double X1 = stats0->GetX1NDC();
      double Y1 = stats0->GetY1NDC();
      double X2 = stats0->GetX2NDC();
      double Y2 = stats0->GetY2NDC();
      WLSEventAction::Instance()->GetHistT(1,SiPM)->SetAxisRange(minbin-5,maxbin+5);
      WLSEventAction::Instance()->GetHistT(1,SiPM)->GetXaxis()->SetTitle("t [ns]");
      WLSEventAction::Instance()->GetHistT(1,SiPM)->SetLineColor(2);
      WLSEventAction::Instance()->GetHistT(1,SiPM)->Draw();
      gPad->Update();
      TPaveStats *stats1 = (TPaveStats*)WLSEventAction::Instance()->GetHistT(1,SiPM)->FindObject("stats");
      stats1->SetTextColor(2);
      stats1->SetLineColor(2);
      stats1->SetX1NDC(X1);
      stats1->SetY1NDC(Y1-(Y2-Y1));
      stats1->SetX2NDC(X2);
      stats1->SetY2NDC(Y1);
      WLSEventAction::Instance()->GetHistT(0,SiPM)->Draw();
      WLSEventAction::Instance()->GetHistT(1,SiPM)->Draw("same");
    }      
    c2.SaveAs("ArrivalTimes.C");
  }
}
