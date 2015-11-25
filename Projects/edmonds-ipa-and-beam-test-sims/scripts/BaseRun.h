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
  static double n_POT_per_microbunch;
  static double n_stopped_muons_per_POT;
  static double n_captured_muons_per_stopped_muon;
  static double n_decayed_muons_per_stopped_muon;  

  static double gas_ionisation_energy;
  static double charge_per_electron;
  static double gas_gain;

 private:
  double fChargeDepositPerHit;

  TH3F* fHitMap;
  TH1F* fSHEDepPlot;
  TH3F* fHitMapTime;

  bool fExists;
};

double BaseRun::n_POT_per_microbunch = 3.1e7;
double BaseRun::n_stopped_muons_per_POT = 0.0019;	   
double BaseRun::n_captured_muons_per_stopped_muon = 0.609;
double BaseRun::n_decayed_muons_per_stopped_muon = 0.391; 
double BaseRun::gas_ionisation_energy = 27.58e-6; // MeV
double BaseRun::charge_per_electron = 1.6e-19;
double BaseRun::gas_gain = 4e4; //  from gas gain fit

BaseRun::BaseRun(std::string filename, std::string particle) : fParticleName(particle), fPDGId(0), 
  fHitMap(NULL), fSHEDepPlot(NULL), fHitMapTime(NULL), fChargeDepositPerHit(0), fNParticlesPerMicrobunch(0), 
  fExists(true) {

  if (fParticleName == "eMinus" || fParticleName == "Flash" || fParticleName == "DIO") {
    fPDGId = 11;
  }
  else if (fParticleName == "pPlus") {
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

  fNSimulatedParticles = fSHDiagChain->GetEntries();
  std::stringstream cutcmd;
  cutcmd << "mcinfo._pdg" << "==" << fPDGId;
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
    std::string drawcmd = "shid._straw:6*shid._device+shid._sector:mcinfo._len>>"+histname;
    std::stringstream cutcmd; cutcmd << "mcinfo._pdg==" << fPDGId;

    fHitMap = new TH3F(histname.c_str(), "", n_bins,min_len_pos,max_len_pos, n_devices*n_sectors,0,n_devices*n_sectors, n_straws,0,n_straws);
    fHitMap->SetXTitle("Position Along Straw [mm]");
    fHitMap->SetYTitle("6*device number + sector number");
    fHitMap->SetZTitle("straw number");
    
    fSHDiagChain->Draw(drawcmd.c_str(), cutcmd.str().c_str(), "goff");
    fHitMap->Sumw2();
    fHitMap->Scale(1/GetNMicrobunchesSimulated());
  }
  return fHitMap;
}

TH1F* BaseRun::GetSHEDepPlot() {
  if (!fSHEDepPlot) {
    double min_edep = 0;
    double max_edep = 0.2;
    double edep_width = 0.001;
    int n_edep_bins = (max_edep - min_edep) / edep_width;
    std::string histname = "SHEDep_" + fParticleName;
    fSHEDepPlot = new TH1F(histname.c_str(), "", n_edep_bins,min_edep,max_edep);

    std::string drawcmd = "edep>>" + histname;
    std::stringstream cutcmd; cutcmd << "mcinfo._pdg==" << fPDGId;
    fSHDiagChain->Draw(drawcmd.c_str(), cutcmd.str().c_str(), "goff");
  }

  return fSHEDepPlot;
}

TH3F* BaseRun::GetHitMapTime() {

  if (!fHitMapTime) {
    double min_time = 0;
    double max_time = 2000;
    double time_width = 1;
    int n_time_bins = (max_time - min_time) / time_width;

    std::string histname = "HitMapTime_" + fParticleName;
    fHitMapTime = new TH3F(histname.c_str(), "", n_time_bins,min_time,max_time, n_devices*n_sectors,0,n_devices*n_sectors, n_straws,0,n_straws);
    fHitMapTime->SetXTitle("Time [ns]");
    fHitMapTime->SetYTitle("6*device number + sector number");
    fHitMapTime->SetZTitle("straw number");
    
    
    std::string drawcmd = "shid._straw:6*shid._device+shid._sector:time>>" + histname;
    std::stringstream cutcmd; cutcmd << "time < " << max_time << " && mcinfo._pdg==" << fPDGId;
    fSHDiagChain->Draw(drawcmd.c_str(), cutcmd.str().c_str(), "goff");    
    fHitMapTime->Scale(1/GetNMicrobunchesSimulated());
  }

  return fHitMapTime;
}

#endif
