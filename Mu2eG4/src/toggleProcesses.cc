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
#include "Geant4/G4ParticleTable.hh"
#include "Geant4/G4VRestProcess.hh"
#include "Geant4/G4ProcessManager.hh"

#include "Geant4/G4EmParameters.hh"

#include "Geant4/G4MuonMinus.hh"
#include "Geant4/G4MuonMinusCapture.hh"
#include "Geant4/G4MuMinusCapturePrecompound.hh"

#include "Geant4/G4PhysicsListHelper.hh"
#include "Geant4/G4GammaConversionToMuons.hh"
#include "Geant4/G4AnnihiToMuPair.hh"
#include "Geant4/G4eeToHadrons.hh"

// Framework includes
#include "cetlib_except/exception.h"

// Mu2e includes
#include "Offline/Mu2eG4/inc/toggleProcesses.hh"
#include "Offline/Mu2eG4/inc/Mu2eRecorderProcess.hh"
#include "Offline/Mu2eG4/inc/Mu2eSpecialCutsProcess.hh"

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
  void addUserProcesses(const Mu2eG4Config::Physics& phys,
                        const Mu2eG4Config::Debug& debug,
                        const Mu2eG4ResourceLimits& lim
                        ) {
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
    // and special cuts process
    if (diagLevel>0) {
      G4cout << __func__
             << " adding mu2eRecorderProcess and "
             << " adding mu2eSpecialCutsProcess to all G4ParticleTable particles" << G4endl;
    }

    Mu2eRecorderProcess* rmp = new Mu2eRecorderProcess();
    rmp->SetVerboseLevel(debug.steppingVerbosityLevel()-1);
    Mu2eSpecialCutsProcess* scp = new Mu2eSpecialCutsProcess(lim);
    scp->SetVerboseLevel(debug.steppingVerbosityLevel()-1);

    G4ParticleTable* ptable = G4ParticleTable::GetParticleTable();
    G4ParticleTable::G4PTblDicIterator* iter = ptable->GetIterator();
    iter->reset();
    while( (*iter)() ){
      G4ParticleDefinition* particle = iter->value();
      G4ProcessManager* pmanager     = particle->GetProcessManager();
      // The process manager takes ownership of the process
      // ph->RegisterProcess(rmp, proton); // RegisterProcess only works for known geant4 processes

      // inline G4int AddRestProcess(G4VProcess* aProcess, G4int ord=ordDefault);
      // inline G4int AddDiscreteProcess(G4VProcess* aProcess, G4int ord=ordDefault);
      // inline G4int AddContinuousProcess(G4VProcess* aProcess, G4int ord=ordDefault);

      // void SetProcessOrderingToLast(
      //                               G4VProcess* aProcess,
      //                               G4ProcessVectorDoItIndex idDoIt
      //                               );
      // Set ordering parameter to the last of all processes

      // enum G4ProcessVectorDoItIndex
      // {
      //   idxAll = -1,         // for all DoIt/GPIL
      //   idxAtRest = 0,       // for AtRestDoIt/GPIL
      //   idxAlongStep = 1,    // for AlongStepDoIt/GPIL
      //   idxPostStep =2,      // for AlongSTepDoIt/GPIL
      //   NDoit =3
      // };

      //G4int AddProcess( G4VProcess* aProcess,
      //                  G4int ordAtRestDoIt = ordInActive,
      //                  G4int ordAlongSteptDoIt = ordInActive,
      //                  G4int ordPostStepDoIt = ordInActive );
      // Adds a process to the process List
      // Return values is the index to the List. Negative return value
      // indicates that the process has not been added due to some errors
      // The first argument is a pointer to the process.
      // Successive arguments are ordering parameters of the process in
      // process vectors. If value is negative, the process is
      // not added to the corresponding process vector

      // enum G4ProcessVectorOrdering
      // {
      //   ordInActive = -1,    // ordering parameter to indicate InActive DoIt
      //   ordDefault = 1000,   // default ordering parameter
      //   ordLast    = 9999    // ordering parameter to indicate the last DoIt
      // };

      // findoing out if the particle can decay (e.g. radioactively or
      // otherwise) to avoid forcing a single process;
      // perhaps checking if there is a AtRestProcess already would be
      // more general

      G4ProcessVector * pVector  = pmanager->GetProcessList();
      G4VProcess *adecayProcess = nullptr;
      for( G4int j=0; j<pmanager->GetProcessListLength(); ++j ) {
        if( ((*pVector)[j]->GetProcessName()).find("Decay") != std::string::npos ) {
          adecayProcess = (*pVector)[j];
          break;
        }
      }
      if (adecayProcess == nullptr || particle->GetPDGLifeTime() <0.0 || particle->GetPDGMass() <= 0.0) {
        // particle is stable
        pmanager->AddProcess(scp, ordInActive, ordInActive, ordDefault);
      } else {
        adecayProcess->SetVerboseLevel(debug.steppingVerbosityLevel()-1);
        pmanager->AddProcess(scp, ordDefault, ordInActive, ordDefault);
      }

      pmanager->AddContinuousProcess(rmp);
      pmanager->SetProcessOrderingToLast(rmp, idxAlongStep); // enum G4ProcessVectorDoItIndex

    }
  }

  //================================================================
} // end namespace mu2e
