#ifndef RANDOMLIMITEDEXPO_HH
#define RANDOMLIMITEDEXPO_HH

//  Extend RandomExponential binding the output variable in a user defined range
// 
// Original author Gianni Onorato
//
//

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandomEngine.h"
#include "FWCore/Services/interface/RandomNumberGeneratorService.h"

namespace mu2e { 
  
  class RandomLimitedExpo {

  public:
    
    RandomLimitedExpo(edm::RandomNumberGeneratorService::base_engine_t& engine);
    
    ~RandomLimitedExpo();
    
    double fire(double tmin, double tmax, double tau = 1);
    
  private:

    CLHEP:: RandFlat _randFlat;


  };

}

#endif
