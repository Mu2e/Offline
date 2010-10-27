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
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"


//CLHEP includes
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Random/RandExponential.h"

namespace mu2e {

  class FoilParticleGenerator {

    friend class ConversionGun;
    friend class DecayInOrbitGun;
    friend class EjectedProtonGun;

  public: 
    
    FoilParticleGenerator( edm::RandomNumberGeneratorService::base_engine_t& engine );
    
    ~FoilParticleGenerator();

    void generatePositionAndTime(CLHEP::Hep3Vector& pos, double& time);

  private:

    // Random numbers generators
    CLHEP::RandFlat     _randFlat;
    CLHEP::RandExponential     _randTime;
    CLHEP::RandGeneral  _randFoils;

    //number of foils
    int _nfoils;

    //Spatial and temporal variables
    double _FPGczmin, _FPGczmax, _FPGphimin, _FPGphimax, _FPGtmin, _FPGtmax;

    //Range for the above
    double _dcz, _dphi, _dt;

    //Build a binned representation of foils volume
    std::vector<double> binnedFoilsVolume();


  };
} // end namespace mu2e,

#endif

