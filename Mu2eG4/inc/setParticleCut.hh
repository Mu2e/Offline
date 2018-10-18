#ifndef Mu2eG4_setParticleCut_hh
#define Mu2eG4_setParticleCut_hh
//
// Set the G4 minimum range/production cut as specified in the config file
//
// Author: KLG based on setMinimumRangeCut
//
//-----------------------------------------------------------------------------

namespace fhicl { class ParameterSet; }

class G4VUserPhysicsList;

namespace mu2e{

  class SimpleConfig;

  void setParticleCut(SimpleConfig const& config, 
                       std::string const& pName, G4VUserPhysicsList* mPL );
  void setParticleCut(fhicl::ParameterSet const& pset, 
                      std::string const& pName,  G4VUserPhysicsList* mPL);

}  // end namespace mu2e

#endif /* Mu2eG4_setParticleCut_hh */
