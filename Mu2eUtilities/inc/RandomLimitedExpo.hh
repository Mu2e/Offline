#ifndef Mu2eUtilities_RandomLimitedExpo_hh
#define Mu2eUtilities_RandomLimitedExpo_hh

//  Extend RandomExponential binding the output variable in a user defined range
//
// Original author Gianni Onorato
//
//

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandomEngine.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

namespace mu2e {

  class RandomLimitedExpo {

  public:

    RandomLimitedExpo(art::RandomNumberGenerator::base_engine_t& engine);

    ~RandomLimitedExpo();

    double fire(double tmin, double tmax, double tau = 1);

  private:

    CLHEP:: RandFlat _randFlat;


  };

}

#endif /* Mu2eUtilities_RandomLimitedExpo_hh */
