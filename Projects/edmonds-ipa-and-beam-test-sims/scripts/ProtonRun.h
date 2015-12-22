#ifndef ProtonRun_h_
#define ProtonRun_h_

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

class ProtonRun : public BaseRun {
 public:
  ProtonRun(std::string filename);

 public:
  //  void Draw();
  //  void DrawProtonsThatReachTracker();
};

ProtonRun::ProtonRun(std::string filename) : BaseRun(filename, "Proton") {

  fChargeDepositPerHit = 8.8e-12 * (BaseRun::gas_gain/4e4); // 4e4 is the gas gain we measured at 1375 V (NB this relationship doesn't quite match the data we have at the moment)
  double n_protons_per_captured_muon = 0.05;
  SetNParticlesPerMicrobunch(n_POT_per_microbunch * n_stopped_muons_per_POT * n_captured_muons_per_stopped_muon * n_protons_per_captured_muon);
}

/*void ProtonRun::DrawProtonsThatReachTracker() {

  double min_energy = 0;
  double max_energy = 50;
  double energy_width = 0.1;
  int n_energy_bins = (max_energy - min_energy) / energy_width;

  TH1F* hKineticEnergy_Gen = new TH1F("hKineticEnergy_Gen", "35M Initial Protons (TDR Cone IPA)", n_energy_bins,min_energy,max_energy);
  GetTrkDiagChain()->Draw("(sqrt(mcgen.mom*mcgen.mom + 938*938) - 938)>>hKineticEnergy_Gen", "mc.pdg==2212", "goff");
  hKineticEnergy_Gen->SetLineWidth(2);
  hKineticEnergy_Gen->SetXTitle("Kinetic Energy [MeV]");
  hKineticEnergy_Gen->SetYTitle("Counts per Simulated Proton");
  hKineticEnergy_Gen->SetStats(false);
  hKineticEnergy_Gen->Scale(1/fNSimulatedParticles);
  hKineticEnergy_Gen->GetYaxis()->SetRangeUser(1e-7, 1e-1);

  // Generated proton kinetic energy of protons that reach tracker
  TH1F* hKineticEnergy_Gen_TrkEnt = new TH1F("hKineticEnergy_Gen_TrkEnt", "", n_energy_bins,min_energy,max_energy);
  GetTrkDiagChain()->Draw("(sqrt(mcgen.mom*mcgen.mom + 938*938) - 938)>>hKineticEnergy_Gen_TrkEnt", "mc.pdg==2212 && mc.ndigi!=0", "goff");
  hKineticEnergy_Gen_TrkEnt->SetLineColor(kRed);
  hKineticEnergy_Gen_TrkEnt->SetLineWidth(2);
  hKineticEnergy_Gen_TrkEnt->SetXTitle("Kinetic Energy [MeV]");
  hKineticEnergy_Gen_TrkEnt->SetYTitle("Counts per Simulated Proton");
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
  }*/
#endif
