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
// $Id: toggleProcesses.cc,v 1.4 2011/05/18 02:27:18 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:18 $
//
//-----------------------------------------------------------------------------


// C++ includes.
#include <vector>
#include <cmath>

// Geant4 includes
#include "G4ParticleTable.hh"
#include "G4VRestProcess.hh"
#include "G4ProcessManager.hh"

// Mu2e includes
#include "Mu2eG4/inc/MuonMinusConversionAtRest.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"

using namespace std;

namespace mu2e{


  void switchDecayOff(const SimpleConfig& config) {

    // Read list of particles for which the decay should be switched off
    vector<int> plist;
    if( ! config.hasName("g4.noDecay") ) return;
    config.getVectorInt("g4.noDecay",plist);

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


 void addUserProcesses(const SimpleConfig& config) {

    G4ParticleTable *theParticleTable = G4ParticleTable::GetParticleTable();

    // so far, only muon is defined as a particle for processes to be added
    G4ParticleDefinition* particle = theParticleTable->FindParticle(13);
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4ProcessVector const* pVector      = pmanager->GetProcessList();

    G4VRestProcess *muProcess = 0;
    if ( config.get<bool>( "g4.doMuMinusConversionAtRest", 1) ) {

        // add muMinusConversionAtRest
        muProcess = new muMinusConversionAtRest( config );
        pmanager->AddRestProcess(muProcess);

         // remove muMinusCaptureAtRest
         for( G4int j=0; j<pmanager->GetProcessListLength(); j++ ) {
            G4VProcess* proc = (*pVector)[j];
            G4String name  = proc->GetProcessName();
            if( name == "muMinusCaptureAtRest" ) pmanager->RemoveProcess(proc);
     }

   }

 }



} // end namespace mu2e
