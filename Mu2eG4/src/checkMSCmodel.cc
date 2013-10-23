//
// In case of use of the ITracker, check if the MSC model matches the ITracker requirements
//
// $Id: checkMSCmodel.cc,v 1.3 2013/10/23 17:22:32 genser Exp $
// $Author: genser $
// $Date: 2013/10/23 17:22:32 $
//

#include "Mu2eG4/inc/checkMSCmodel.hh"
#include "ConfigTools/inc/SimpleConfig.hh"

#include "messagefacility/MessageLogger/MessageLogger.h"

// Geant4 includes
#include "G4String.hh"
#include "G4ParticleTable.hh"
#include "G4Run.hh"
#include "G4Timer.hh"
#include "G4ProcessManager.hh"

#include "G4eMultipleScattering.hh"
#if G4VERSION<4095
  #include "G4UrbanMscModel92.hh"
#endif

#include <iostream>

namespace mu2e{

  void checkMSCmodel( SimpleConfig const& config ){

          if (config.getBool("hasITracker",false)) {
          #if G4VERSION<4095
                  bool change = false;
                  if (  config.getBool("itracker.changeMSC",false) ) {
                          G4ParticleTable *theParticleTable = G4ParticleTable::GetParticleTable();
                          G4ParticleDefinition* particle = theParticleTable->FindParticle(11);
                          G4ProcessManager* pmanager = particle->GetProcessManager();
                          G4ProcessVector const* pVector = pmanager->GetProcessList();
                          for( G4int j=0; j<pmanager->GetProcessListLength(); j++ ) {
                                  G4VProcess* proc = (*pVector)[j];
                                  G4String name  = proc->GetProcessName();
                                  if( name == "msc" ) {
                                          pmanager->RemoveProcess(proc);
                                          change = true;
                                  }
                          }
                          if (change) {
                                  G4eMultipleScattering* msc = new G4eMultipleScattering();
                                  msc->AddEmModel(0, new G4UrbanMscModel92());
                                  pmanager->AddProcess(msc,                       -1, 1, 1);
                          }
                  }
          #endif
          }

  }

}  // end namespace mu2e
