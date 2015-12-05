#ifndef DeuteronRun_h_
#define DeuteronRun_h_


#include "BaseRun.h"

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH2.h"
#include "TH3.h"
#include "TLatex.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TLegend.h"

#include <iostream>
#include <sstream>
#include <algorithm>

class DeuteronRun : public BaseRun {
 public:
  DeuteronRun(std::string filename);

 private:
  double fChargeDepositPerHit;

 public:
  //  void Draw();
  void DrawDeuteronsThatReachTracker();
};

DeuteronRun::DeuteronRun(std::string filename) : BaseRun(filename, "Deuteron") {

  //  fChargeDepositPerHit = 8.8e-12;
  //  double n_deuterons_per_captured_muon = 0.0175;
  double n_deuterons_per_captured_muon = 0.025; // CD3 number
  fNParticlesPerMicrobunch = n_POT_per_microbunch * n_stopped_muons_per_POT * n_captured_muons_per_stopped_muon * n_deuterons_per_captured_muon;
}
/*
void DeuteronRun::Draw() {

  int n_straws = 48;
  int n_layers = 2;
  int n_sectors = 6;
  int n_devices = 44;

  int min_full_straw_number = 0;
  int max_full_straw_number = 100;
  int min_device_number = 0;
  int max_device_number = 500;
  double bin_width = 1;
  int n_full_straw_bins = (max_full_straw_number - min_full_straw_number)/bin_width;
  int n_device_bins = (max_device_number - min_device_number)/bin_width;

  TH2F* hHitLocations = new TH2F("hHitLocations", "Straw Hits per Initial Deuteron", n_device_bins,min_device_number,max_device_number, n_full_straw_bins,min_full_straw_number,max_full_straw_number);
  hHitLocations->SetXTitle("Device*10 + Sector");
  hHitLocations->SetYTitle("Straw Number");
  hHitLocations->GetZaxis()->SetLabelSize(0.03);
  hHitLocations->SetStats(false);
  //  hHitLocations->SetZTitle("Hits per Emitted Deuteron");

  fSHDiagChain->Draw("(straw*2 + layer):(device*10 + sector)>>hHitLocations", "mcpdg==1000010020 && time>700 && time<1695", "COLZ"); 
  hHitLocations->Scale(1.0/fNSimulatedParticles);

  std::stringstream latex;
  latex.precision(3);
  latex << "#splitline{Total Number of Straw}{Hits Per Deuteron = " << GetNStrawHitsPerSimulatedParticle() << " #pm ";
  latex.precision(1);
  latex << GetNStrawHitsPerSimulatedParticleError() << "}";
  
  TLatex* text = new TLatex(1, 80, latex.str().c_str());
  text->SetTextSize(0.05);
  text->SetTextColor(kRed);
  text->Draw("SAME");
}
*/
void DeuteronRun::DrawDeuteronsThatReachTracker() {

  double min_energy = 0;
  double max_energy = 50;
  double energy_width = 0.1;
  int n_energy_bins = (max_energy - min_energy) / energy_width;

  TH1F* hKineticEnergy_Gen = new TH1F("hKineticEnergy_Gen", "35M Initial Deuterons (TDR Cone IPA)", n_energy_bins,min_energy,max_energy);
  fTrkDiagChain->Draw("(sqrt(mcgen.mom*mcgen.mom + 938*938) - 938)>>hKineticEnergy_Gen", "mc.pdg==1000010020", "goff");
  hKineticEnergy_Gen->SetLineWidth(2);
  hKineticEnergy_Gen->SetXTitle("Kinetic Energy [MeV]");
  hKineticEnergy_Gen->SetYTitle("Counts per Simulated Deuteron");
  hKineticEnergy_Gen->SetStats(false);
  hKineticEnergy_Gen->Scale(1/fNSimulatedParticles);
  hKineticEnergy_Gen->GetYaxis()->SetRangeUser(1e-7, 1e-1);

  // Generated proton kinetic energy of protons that reach tracker
  TH1F* hKineticEnergy_Gen_TrkEnt = new TH1F("hKineticEnergy_Gen_TrkEnt", "", n_energy_bins,min_energy,max_energy);
  fTrkDiagChain->Draw("(sqrt(mcgen.mom*mcgen.mom + 938*938) - 938)>>hKineticEnergy_Gen_TrkEnt", "mc.pdg==1000010020 && mc.ndigi!=0", "goff");
  hKineticEnergy_Gen_TrkEnt->SetLineColor(kRed);
  hKineticEnergy_Gen_TrkEnt->SetLineWidth(2);
  hKineticEnergy_Gen_TrkEnt->SetXTitle("Kinetic Energy [MeV]");
  hKineticEnergy_Gen_TrkEnt->SetYTitle("Counts per Simulated Deuteron");
  hKineticEnergy_Gen_TrkEnt->SetStats(false);
  hKineticEnergy_Gen_TrkEnt->Scale(1/fNSimulatedParticles);
  hKineticEnergy_Gen_TrkEnt->GetYaxis()->SetRangeUser(1e-7, 1e-1);


  TCanvas* c1 = new TCanvas("c1", "c1");
  c1->SetLogy();
  hKineticEnergy_Gen->Draw();
  hKineticEnergy_Gen_TrkEnt->Draw("SAME");

  TLegend* legend = new TLegend(0.42, 0.85, 0.72, 0.65, "");
  legend->SetBorderSize(0);
  legend->SetTextSize(0.035);
  legend->SetFillColor(kWhite);
  legend->AddEntry(hKineticEnergy_Gen, "Kinetic Energy at Emission (all)", "l");
  legend->AddEntry(hKineticEnergy_Gen_TrkEnt, "Kinetic Energy at Emission (reach tracker)", "l");
  legend->Draw();
}
#endif
