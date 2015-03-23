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
} // end namespace mu2e
