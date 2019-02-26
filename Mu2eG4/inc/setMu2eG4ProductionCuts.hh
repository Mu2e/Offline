#ifndef Mu2eG4_setMu2eG4ProductionCuts_hh
#define Mu2eG4_setMu2eG4ProductionCuts_hh
//
// Set the G4 range/production cuts including the ones for regions if any
// Original author  K.L.Genser

namespace fhicl { class ParameterSet; }

class G4VModularPhysicsList;

namespace mu2e{

  void setMu2eG4ProductionCuts(const fhicl::ParameterSet& pset, G4VModularPhysicsList* mPL);

}  // end namespace mu2e

#endif /* Mu2eG4_setMu2eG4ProductionCuts_hh */
