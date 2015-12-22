#ifndef BaseRun_h_
#define BaseRun_h_
// A base class for all useful runs for different types of particles

#include <iostream>
#include <sstream>
#include <stdexcept>

#include "TFile.h"
#include "TChain.h"

#include "TH1.h"
#include "TH3.h"

class BaseRun {

 public:
  BaseRun(std::string filename, std::string particle);

  std::string GetParticleName() { return fParticleName; }
  double GetNSimulatedParticles() { return fNSimulatedParticles; }
  double GetNParticlesPerMicrobunch() { return fNParticlesPerMicrobunch; }
  double GetNMicrobunchesSimulated() { return fNSimulatedParticles/fNParticlesPerMicrobunch; }
  double GetNStrawHits() { return fNStrawHits; }
  double GetNStrawHitsError() { return fNStrawHitsError; }
  double GetChargeDepositPerHit() { return fChargeDepositPerHit; }

  double GetNStrawHitsPerSimulatedParticle() { return fNStrawHits/fNSimulatedParticles; }
  double GetNStrawHitsPerSimulatedParticleError() { return fNStrawHitsError/fNSimulatedParticles; }

  TH1F* GetSHEDepPlot();
  TH3F* GetHitMap();

  TH3F* GetHitMapTime();
  TH3F* GetHitMapEDep();

  bool Exists() { return fExists; }

  TChain* GetTrkDiagChain() { return fTrkDiagChain; }
  TChain* GetSHDiagChain() { return fSHDiagChain; }

  void SetNParticlesPerMicrobunch(double n_particles) { fNParticlesPerMicrobunch = n_particles; }

 private:
  std::string fParticleName;
  int fPDGId;

  TChain* fTrkDiagChain;
  TChain* fSHDiagChain;

  double fNSimulatedParticles;
  double fNParticlesPerMicrobunch;
  double fNStrawHits;
  double fNStrawHitsError;

  static const int n_devices = 44;
  static const int n_sectors = 6;
  //  static const int n_layers = 2;
  static const int n_straws = 96;

 public:
  static double n_years_running;
  static double n_POT_per_year;
  static double n_POT_per_microbunch;
  static double n_microbunches_per_year;
  static double n_stopped_muons_per_POT;
  static double n_captured_muons_per_stopped_muon;
  static double n_decayed_muons_per_stopped_muon;  

  static double gas_ionisation_energy;
  static double charge_per_electron;
  static double gas_gain;

  static double max_microbunch_time;

 private:
  double fChargeDepositPerHit;

  TH3F* fHitMap;
  TH1F* fSHEDepPlot;
  TH3F* fHitMapTime;
  TH3F* fHitMapEDep;

  bool fExists;
};

double BaseRun::max_microbunch_time = 1700;

double BaseRun::n_years_running = 3;
double BaseRun::n_POT_per_year = 1.2e20;
//double BaseRun::n_POT_per_microbunch = 3.1e7; //TDR
double BaseRun::n_POT_per_microbunch = 3.9e7; //CD3
double BaseRun::n_microbunches_per_year = BaseRun::n_POT_per_year / BaseRun::n_POT_per_microbunch;
double BaseRun::n_stopped_muons_per_POT = 0.0019;	   
double BaseRun::n_captured_muons_per_stopped_muon = 0.609;
double BaseRun::n_decayed_muons_per_stopped_muon = 0.391; 
double BaseRun::gas_ionisation_energy = 27.58e-6; // MeV
double BaseRun::charge_per_electron = 1.6e-19;
double BaseRun::gas_gain = 4e4; //  from gas gain fit
//double BaseRun::gas_gain = 7e4; // new CD3 possibility

BaseRun::BaseRun(std::string filename, std::string particle) : fParticleName(particle), fPDGId(0), 
  fHitMap(NULL), fSHEDepPlot(NULL), fHitMapTime(NULL), fHitMapEDep(NULL), fChargeDepositPerHit(0), fNParticlesPerMicrobunch(0), 
  fExists(true) {

  if (fParticleName == "ConvEMinus" || fParticleName == "Flash" || fParticleName == "DIO") {
    fPDGId = 11;
  }
  else if (fParticleName == "Proton") {
    fPDGId = 2212;
  }
  else if (fParticleName == "Deuteron") {
    fPDGId = 1000010020;
  }

  if (fPDGId == 0) {
    std::string message = "ERROR: BaseRun::BaseRun() : Particle " + particle + " is not suppored in BaseRun";
    throw std::invalid_argument(message);
  }

  // Get the shdiag tree and chain multiple files together
  std::string trkdiag_name = "TRFDownstreameMinus/trkdiag"; // all runs are called TRFDownstreameMinus
  fTrkDiagChain = new TChain(trkdiag_name.c_str());
  fTrkDiagChain->Add(filename.c_str());
  if (fTrkDiagChain->GetEntries() == 0) {
    std::string message = "ERROR: BaseRun::BaseRun() : TrkDiagTree " + trkdiag_name + " has no entries (does it exist in this file?)";
    //    throw std::invalid_argument(message);
    fExists = false;
  }

  // Get the shdiag tree and chain multiple files together
  std::string shdiag_name = "makeSH/shdiag";
  fSHDiagChain = new TChain(shdiag_name.c_str());
  fSHDiagChain->Add(filename.c_str());
  if (fSHDiagChain->GetEntries() == 0) {
    std::string message = "ERROR: BaseRun::BaseRun() : ShDiagTree " + shdiag_name + " has no entries (does it exist in this file?)";
    //    throw std::invalid_argument(message);
    fExists = false;
  }

  fNSimulatedParticles = fTrkDiagChain->GetEntries();
  std::stringstream cutcmd; cutcmd << ""; //cutcmd << "mcinfo._pdg==" << fPDGId;
  fNStrawHits = fSHDiagChain->GetEntries(cutcmd.str().c_str());
  fNStrawHitsError = std::sqrt(fNStrawHits);

  fChargeDepositPerHit = (GetSHEDepPlot()->GetMean() / gas_ionisation_energy) * charge_per_electron * gas_gain;
  //  std::cout << filename << " " << fParticleName << " " << fSHDiagChain->GetEntries() << " " << fSHDiagChain->GetEntries() << " worked!" << std::endl;		   
}

TH3F* BaseRun::GetHitMap() {

  if (!fHitMap) {
    // Create the hit map
    double min_len_pos = -600;
    double max_len_pos = 600;
    double bin_width = 10; // mm
    int n_bins = (max_len_pos - min_len_pos) / bin_width;
    
    std::string histname = "fHitMap_" + fParticleName;
    std::string drawcmd;
    if (fParticleName != "ConvEMinus") {
      drawcmd = "shid._straw:6*shid._device+shid._sector:mcinfo._len>>"+histname;
    }
    else {
      drawcmd = "tsh._straw:6*tsh._device+tsh._sector:tshmc._len>>"+histname;
    }

    fHitMap = new TH3F(histname.c_str(), "", n_bins,min_len_pos,max_len_pos, n_devices*n_sectors,0,n_devices*n_sectors, n_straws,0,n_straws);
    fHitMap->SetXTitle("Position Along Straw [mm]");
    fHitMap->SetYTitle("6*device number + sector number");
    fHitMap->SetZTitle("straw number");
    
    if (fParticleName != "ConvEMinus") {
      fSHDiagChain->Draw(drawcmd.c_str(), "", "goff");
    }
    else {
      fTrkDiagChain->Draw(drawcmd.c_str(), "", "goff");
    }
    fHitMap->Sumw2();
    fHitMap->Scale(1/GetNMicrobunchesSimulated());
  }
  return fHitMap;
}

TH3F* BaseRun::GetHitMapEDep() {

  if (!fHitMapEDep) {
    // Create the hit map
    double min_edep = 0;
    double max_edep = 0.02;
    double bin_width = 0.0001; // MeV
    int n_bins = (max_edep - min_edep) / bin_width;
    
    std::string histname = "fHitMapEDep_" + fParticleName;
    std::string drawcmd;
    if (fParticleName != "ConvEMinus") {
      drawcmd = "shid._straw:6*shid._device+shid._sector:edep>>"+histname;
    }
    else {
      drawcmd = "tsh._straw:6*tsh._device+tsh._sector:tsh._edep>>"+histname;
    }

    fHitMapEDep = new TH3F(histname.c_str(), "", n_bins,min_edep,max_edep, n_devices*n_sectors,0,n_devices*n_sectors, n_straws,0,n_straws);
    fHitMapEDep->SetXTitle("Energy Depositied [MeV]");
    fHitMapEDep->SetYTitle("6*device number + sector number");
    fHitMapEDep->SetZTitle("straw number");
    
    if (fParticleName != "ConvEMinus") {
      fSHDiagChain->Draw(drawcmd.c_str(), "", "goff");
    }
    else {
      fTrkDiagChain->Draw(drawcmd.c_str(), "", "goff");
    }
    fHitMapEDep->Sumw2();
    fHitMapEDep->Scale(1/GetNMicrobunchesSimulated());
  }
  return fHitMapEDep;
}

TH1F* BaseRun::GetSHEDepPlot() {
  if (!fSHEDepPlot) {
    double min_edep = 0;
    double max_edep = 0.2;
    double edep_width = 0.001;
    int n_edep_bins = (max_edep - min_edep) / edep_width;
    std::string histname = "SHEDep_" + fParticleName + "Run";
    fSHEDepPlot = new TH1F(histname.c_str(), "Energy Deposited by All Particles", n_edep_bins,min_edep,max_edep);

    std::string drawcmd = "mcinfo._energy>>" + histname;
    fSHDiagChain->Draw(drawcmd.c_str(), "", "goff");
  }

  return fSHEDepPlot;
}

TH3F* BaseRun::GetHitMapTime() {

  if (!fHitMapTime) {
    double min_time = 0;
    double max_time = BaseRun::max_microbunch_time;
    double time_width = 1;
    int n_time_bins = (max_time - min_time) / time_width;

    std::string histname = "HitMapTime_" + fParticleName;
    fHitMapTime = new TH3F(histname.c_str(), "", n_time_bins,min_time,max_time, n_devices*n_sectors,0,n_devices*n_sectors, n_straws,0,n_straws);
    fHitMapTime->SetXTitle("Time [ns]");
    fHitMapTime->SetYTitle("6*device number + sector number");
    fHitMapTime->SetZTitle("straw number");
    
    
    std::string drawcmd = "shid._straw:6*shid._device+shid._sector:time>>" + histname;
    std::stringstream cutcmd; cutcmd << "time < " << max_time;
    fSHDiagChain->Draw(drawcmd.c_str(), cutcmd.str().c_str(), "goff");    
    fHitMapTime->Scale(1/GetNMicrobunchesSimulated());
  }

  return fHitMapTime;
}

#endif
