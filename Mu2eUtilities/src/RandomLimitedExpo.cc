//
// Extend RandomExponential binding the output variable in a user defined range
//
// Original author Gianni Onorato
//

//C++ includes
#include <cmath>

//CLHEP includes
#include "CLHEP/Random/RandFlat.h"


// Framework includes
#include "Mu2eUtilities/inc/RandomLimitedExpo.hh"

namespace mu2e{

  RandomLimitedExpo::RandomLimitedExpo(art::RandomNumberGenerator::base_engine_t& engine):
    _randFlat(engine)
  {}

  RandomLimitedExpo::~RandomLimitedExpo()
  {}

  double RandomLimitedExpo::fire(double tmin, double tmax, double tau) {
    double lambda = 1/tau;
    double flatlimit = 1 - exp(lambda*(tmin-tmax));
    double randvar = (flatlimit*_randFlat.fire());
    return (-log(1-randvar)*tau)+tmin;
  }
}
