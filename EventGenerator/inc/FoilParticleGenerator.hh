#ifndef FOILPARTICLEGENERATOR_HH
#define FOILPARTICLEGENERATOR_HH

//
// Generate position and time of a generic particle coming from target foils.
//
//
// For now this is limited to:
//  - Uniform over the targets.
//  - Need to be improved at a later date.
//  - Limits on cos(theta) and phi but uniform within the range.
 

#include <memory>

// Framework includes
#include "FWCore/Services/interface/RandomNumberGeneratorService.h"

// Mu2e includes
#include "Mu2eUtilities/inc/PDGCode.hh"
#include "Mu2eUtilities/inc/RandomLimitedExpo.hh"

//CLHEP includes
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Random/RandExponential.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class FoilParticleGenerator {

  public: 
    
    FoilParticleGenerator( edm::RandomNumberGeneratorService::base_engine_t& engine,
                           double tmin, double tmax );

    ~FoilParticleGenerator();

    void generatePositionAndTime(CLHEP::Hep3Vector& pos, double& time);

  private:

    // time generation range
    double _tmin, _tmax;

    //number of foils
    int _nfoils;

    // Random numbers generators
    CLHEP::RandFlat     _randFlat;
    //    CLHEP::RandExponential     _randTime;
    RandomLimitedExpo _randTime;
    CLHEP::RandGeneral  _randFoils;

    //Build a binned representation of foils volume
    std::vector<double> binnedFoilsVolume();


  };
} // end namespace mu2e,

#endif

