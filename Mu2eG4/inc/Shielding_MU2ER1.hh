// Physics list inheriting from Shielding using regions etc.
// Uses ShieldingM
// Original author K.L. Genser January 2019
//
#ifndef Mu2eG4_Shielding_MU2ER1_h
#define Mu2eG4_Shielding_MU2ER1_h 1

#include "Shielding.hh"

namespace mu2e {

  class Shielding_MU2ER1: public Shielding // could try to templetize it here
  {
  public:

    explicit Shielding_MU2ER1( const fhicl::ParameterSet& pSet, G4int verbose = 1,
                               G4String low_energy_neutron_model = "HP",
                               G4String HadrPhysVariant = "M");

    virtual ~Shielding_MU2ER1(){};

    virtual void SetCuts() override;

  private:

    const fhicl::ParameterSet& pset;

  };

}
#include "Shielding_MU2ER1.icc"

#endif
