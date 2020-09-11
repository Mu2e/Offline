//
// Class Description:
//
// Function that handles the switching on and off of G4 processes.  This
// is handled through the configuration files and includes the following
// commands:
//
// g4.noDecay - turns off decays of specified particles
// muMinusConversionAtRest.do - turns on the at rest G4 process
// MuonMinusConversionAtRest and turns off MuonMinusCaptureAtRest
//
//
//-----------------------------------------------------------------------------


// C++ includes.
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

// Geant4 includes
#include "G4ParticleTable.hh"
#include "G4VRestProcess.hh"
#include "G4ProcessManager.hh"

#include "G4EmParameters.hh"

#include "G4MuonMinus.hh"
#include "G4MuonMinusCapture.hh"
#include "G4MuMinusCapturePrecompound.hh"

#include "G4PhysicsListHelper.hh"
#include "G4GammaConversionToMuons.hh"
#include "G4AnnihiToMuPair.hh"
#include "G4eeToHadrons.hh"

// Framework includes
#include "cetlib_except/exception.h"

// Mu2e includes
#include "Mu2eG4/inc/Mu2eRecorderProcess.hh"

#include "Mu2eG4/inc/toggleProcesses.hh"

namespace mu2e{

  //================================================================
  void switchDecayOff(const Mu2eG4Config::Physics& phys, const Mu2eG4Config::Debug& debug) {

    std::vector<int> plist = phys.noDecay();
    int diagLevel = debug.diagLevel();

    G4ParticleTable *theParticleTable = G4ParticleTable::GetParticleTable();
    for( size_t i=0; i<plist.size(); ++i ) {
      int pdg = plist[i];
      G4ParticleDefinition* particle = theParticleTable->FindParticle(pdg);
      if( particle==nullptr && diagLevel>-1 ) {
        G4cout << __func__ << " : cannot find particle pdgId=" << pdg << G4endl;
      } else {
        G4ProcessManager* pmanager = particle->GetProcessManager();
        G4ProcessVector * pVector  = pmanager->GetProcessList();
        G4VProcess *decayProcess = nullptr;
        for( G4int j=0; j<pmanager->GetProcessListLength(); ++j ) {
          if( (*pVector)[j]->GetProcessName() == "Decay" ) {
            decayProcess = (*pVector)[j];
            break;
          }
        }
        if( decayProcess==nullptr && diagLevel>-1 ) {
          G4cout << __func__ << " : cannot find decay process for particle pdgId=" << pdg
               << " (" << particle->GetParticleName() << ")" << G4endl;
        } else {
          pmanager->RemoveProcess(decayProcess);
          if (diagLevel>0) {
            G4cout << __func__ << " : decay process is removed for particle pdgId=" << pdg
                   << " (" << particle->GetParticleName() << ")" << G4endl;
          }
        }
        if (diagLevel>0) {
          G4cout << __func__ << " : list of processes defined for particle pdgId=" << pdg
                 << " (" << particle->GetParticleName() << "):" << G4endl;
          for( G4int j=0; j<pmanager->GetProcessListLength(); ++j )
            G4cout << (*pVector)[j]->GetProcessName() << G4endl;
        }
      }
    }

  }

  void switchCaptureDModel(const Mu2eG4Config::Physics& phys, const Mu2eG4Config::Debug& debug) {

    std::string cDModel= phys.captureDModel();
    int diagLevel = debug.diagLevel();

    // change muMinusCaptureAtRest deexcitation model to muMinusNuclearCapture
    // this is limited to muon minus only

    if (std::string("muMinusNuclearCapture")!=cDModel) {
      return;
    }

    G4ParticleDefinition* particle = G4MuonMinus::MuonMinus();
    if( particle==0 && diagLevel>-1 ) {
      G4cout << __func__ << " : cannot find MuonMinus " << G4endl;
    } else {
      G4ProcessManager* pmanager = particle->GetProcessManager();
      G4ProcessVector * pVector  = pmanager->GetProcessList();
      G4VProcess* muCapProcess = nullptr;
      for( G4int j=0; j<pmanager->GetProcessListLength(); ++j ) {
        if( (*pVector)[j]->GetProcessName() == "muMinusCaptureAtRest" ) {
          muCapProcess = (*pVector)[j];
          break;
        }
      }
      if( muCapProcess==0 && diagLevel>-1 ) {
        G4cout << __func__ << " : cannot find muMinusCaptureAtRest process for "
             << particle->GetParticleName() << G4endl;
      } else {
        pmanager->RemoveProcess(muCapProcess);
        if (diagLevel>0) {
          G4cout << __func__ << " : muMinusCaptureAtRest process is removed for "
                 << particle->GetParticleName() << G4endl;
        }
        G4MuonMinusCapture* muProcess = new G4MuonMinusCapture(new G4MuMinusCapturePrecompound());
        G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
        ph->RegisterProcess(muProcess, particle);
        //        pmanager->AddRestProcess( muProcess );
        if (diagLevel>0) {
          G4cout << __func__ << " : added muMinusCaptureAtRest with muMinusNuclearCapture for "
                 << particle->GetParticleName() << G4endl;
        }
      }
      if (diagLevel>0) {
        G4cout  << __func__ << " : list of processes defined for "
                << particle->GetParticleName() << " :" << G4endl;
        for( G4int j=0; j<pmanager->GetProcessListLength(); ++j ) {
          G4cout << (*pVector)[j]->GetProcessName() << G4endl;
        }
      }
    }
  }

  //================================================================
  void addUserProcesses(const Mu2eG4Config::Physics& phys, const Mu2eG4Config::Debug& debug) {
    std::vector<std::string> plist = phys.addProcesses();

    int diagLevel = debug.diagLevel();

    // search for process names, we assume process implies specific particle, for now

    // for ( const auto& elem : plist) {
    //   G4cout << __func__ << " " << elem << G4endl;
    // }

    G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

    G4ParticleDefinition* gamma = G4Gamma::Gamma();
    G4ParticleDefinition* positron = G4Positron::Positron();

    // consider replaceing implementation with one using MCDataProducts/inc/ProcessCode

    std::vector<std::string>::const_iterator ipos;

    ipos = std::find(plist.begin(), plist.end(), std::string("GammaToMuPair"));
    bool gmumu = ( ipos != plist.end() );

    ipos = std::find(plist.begin(), plist.end(), std::string("AnnihiToMuPair"));
    bool pmumu = ( ipos != plist.end() );

    ipos = std::find(plist.begin(), plist.end(), std::string("ee2hadr"));
    bool phad = ( ipos != plist.end() );

    // processes are owned by process managers
    if(gmumu) {
      if (diagLevel>0) {
        G4cout << __func__ << " adding GammaToMuPair process to gamma" << G4endl;
      }
      G4GammaConversionToMuons* theGammaToMuMu = new G4GammaConversionToMuons();
      ph->RegisterProcess(theGammaToMuMu, gamma);
    }
    if(pmumu) {
      if (diagLevel>0) {
        G4cout << __func__ << " adding AnnihiToMuPair process to positron" << G4endl;
      }
      G4AnnihiToMuPair* thePosiToMuMu = new G4AnnihiToMuPair();
      ph->RegisterProcess(thePosiToMuMu, positron);
    }
    if(phad) {
      if (diagLevel>0) {
        G4cout << __func__ << " adding ee2hadr process to positron" << G4endl;
      }
      G4eeToHadrons* thePosiToHadrons = new G4eeToHadrons();
      ph->RegisterProcess(thePosiToHadrons, positron);
    }

    // special process to look at the track before post step interaction
    if (diagLevel>0) {
      G4cout << __func__
             << " adding Mu2eRecorderProcess to all G4ParticleTable particles" << G4endl;
    }

    Mu2eRecorderProcess* rmp = new Mu2eRecorderProcess();
    rmp->SetVerboseLevel(debug.steppingVerbosityLevel()-1);
    G4ParticleTable* ptable = G4ParticleTable::GetParticleTable();
    G4ParticleTable::G4PTblDicIterator* iter = ptable->GetIterator();
    iter->reset();
    while( (*iter)() ){
      G4ParticleDefinition* particle = iter->value();
      G4ProcessManager* pmanager     = particle->GetProcessManager();
      // The process manager takes ownership of the process
      // ph->RegisterProcess(rmp, proton); // RegisterProcess only works for known geant4 processes
      pmanager->AddContinuousProcess(rmp);
      pmanager->SetProcessOrderingToLast(rmp, idxAlongStep);
    }
  }

  //================================================================
} // end namespace mu2e
