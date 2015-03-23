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
// $Id: toggleProcesses.cc,v 1.9 2013/03/26 19:39:33 genser Exp $
// $Author: genser $
// $Date: 2013/03/26 19:39:33 $
//
//-----------------------------------------------------------------------------


// C++ includes.
#include <vector>
#include <cmath>

// Geant4 includes
#include "G4ParticleTable.hh"
#include "G4VRestProcess.hh"
#include "G4ProcessManager.hh"

// Framework includes
#include "cetlib/exception.h"

// Mu2e includes
#include "Mu2eG4/inc/MuonMinusConversionAtRest.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "fhiclcpp/ParameterSet.h"

using namespace std;

namespace mu2e{

  //================================================================
  void switchDecayOff(const std::vector<int>& plist) {

    G4ParticleTable *theParticleTable = G4ParticleTable::GetParticleTable();
    for( size_t i=0; i<plist.size(); ++i ) {
      int pdg = plist[i];
      G4ParticleDefinition* particle = theParticleTable->FindParticle(pdg);
      if( particle==0 ) {
        cout << "SwitchDecayOff: cannot find particle pdgId=" << pdg << endl;
      } else {
        G4ProcessManager* pmanager = particle->GetProcessManager();
        G4ProcessVector * pVector  = pmanager->GetProcessList();
        G4VProcess *decayProcess = 0;
        for( G4int j=0; j<pmanager->GetProcessListLength(); j++ ) {
          if( (*pVector)[j]->GetProcessName() == "Decay" ) {
            decayProcess = (*pVector)[j];
            break;
          }
        }
        if( decayProcess==0 ) {
          cout << "SwitchDecayOff: cannot find decay process for particle pdgId=" << pdg
               << " (" << particle->GetParticleName() << ")" << endl;
        } else {
          pmanager->RemoveProcess(decayProcess);
          cout << "SwitchDecayOff: decay process is removed for particle pdgId=" << pdg
               << " (" << particle->GetParticleName() << ")" << endl;
        }
        cout << "SwitchDecayOff: list of processes defined for particle pdgId=" << pdg
             << " (" << particle->GetParticleName() << "):" << endl;
        for( G4int j=0; j<pmanager->GetProcessListLength(); j++ )
          cout << (*pVector)[j]->GetProcessName() << endl;
      }
    }

  }

  void switchDecayOff(const SimpleConfig& config) {
    // Read list of particles for which the decay should be switched off
    std::vector<int> plist;
    if( ! config.hasName("g4.noDecay") ) return;
    config.getVectorInt("g4.noDecay",plist);
    switchDecayOff(plist);
  }

  void switchDecayOff(const fhicl::ParameterSet& pset) {
    std::vector<int> plist = pset.get<std::vector<int> >("noDecay");
    switchDecayOff(plist);
  }

  //================================================================
  template<class Config>
  void addUserProcesses(bool muConversionAtRest, bool muAtomicCapture, const Config& config) {

    G4ParticleTable *theParticleTable = G4ParticleTable::GetParticleTable();

    // so far, only muon is defined as a particle for processes to be added
    G4ParticleDefinition* particle = theParticleTable->FindParticle(13);
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4ProcessVector const* pVector = pmanager->GetProcessList();
    G4String const & particleName = particle->GetParticleName();

    G4String pNameToRemove("muMinusCaptureAtRest");

    if ( muAtomicCapture && muConversionAtRest )
      throw cet::exception("muon processes")
        << "Can not have both "
        << " g4.useNewMuMinusAtomicCapture and g4.doMuMinusConversionAtRest\n";

    if ( muConversionAtRest ) {

      cout << __func__ << " " << particleName << " processes: " << endl;
      for ( G4int j=0; j<pmanager->GetProcessListLength(); ++j ) {
        G4VProcess* proc = (*pVector)[j];
        G4String const & processName = proc->GetProcessName();
        cout << __func__ << " " << processName << endl;
      }

      // add muMinusConversionAtRest
      G4VRestProcess *muProcess = new muMinusConversionAtRest( config );
      cout << __func__ << " Adding " << muProcess->GetProcessName()
           << " to " << particleName << endl;
      pmanager->AddRestProcess(muProcess);

      // remove muMinusCaptureAtRest
      for ( G4int j=0; j<pmanager->GetProcessListLength(); ++j ) {
        G4VProcess* proc = (*pVector)[j];
        G4String const & processName  = proc->GetProcessName();
        if ( processName == pNameToRemove ) {
          cout << __func__ << " Removing " << processName << " from " << particleName << endl;
          pmanager->RemoveProcess(proc);
          break;
        }
      }
    }
  }

  void addUserProcesses(const SimpleConfig& config) {
    bool muConversionAtRest = config.getBool( "g4.doMuMinusConversionAtRest", false);
    bool muAtomicCapture    = config.getBool( "g4.useNewMuMinusAtomicCapture", false);
    addUserProcesses(muConversionAtRest, muAtomicCapture, config);
  }

  void addUserProcesses(const fhicl::ParameterSet& pset) {
    bool muConversionAtRest = false;
    bool muAtomicCapture    = pset.get<bool>( "useNewMuMinusAtomicCapture");
    addUserProcesses(muConversionAtRest, muAtomicCapture, pset);
  }

  //================================================================
} // end namespace mu2e
