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
#include "TargetGeom/inc/Target.hh"

//CLHEP includes
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Random/RandExponential.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class FoilParticleGenerator {

  public: 
    
    enum foilGen_enum {
      flatFoil, volWeightFoil, expoFoil, expoVolWeightFoil
    };

    enum posGen_enum {
      flatPos
    };

    enum timeGen_enum {
      flatTime, limitedExpoTime
    };



    FoilParticleGenerator( edm::RandomNumberGeneratorService::base_engine_t& engine,
                           double tmin, double tmax, foilGen_enum foilAlgo, 
                           posGen_enum  posAlgo, timeGen_enum  timeAlgo);
    
    ~FoilParticleGenerator();
    
    void generatePositionAndTime(CLHEP::Hep3Vector& pos, double& time); 
    
    int iFoil();

  private:

    // time generation range
    double _tmin, _tmax;

    // foil, position and time random algorithm
    foilGen_enum  _foilAlgo;
    posGen_enum   _posAlgo;
    timeGen_enum  _timeAlgo;

    //number of foils
    int _nfoils;

    //extracted foil.
    int _ifoil;

    // Random numbers generators
    CLHEP::RandFlat     _randFlat;
    RandomLimitedExpo   _randTime;
    CLHEP::RandGeneral  _randFoils;
    CLHEP::RandGeneral  _randExpoFoils;

    //Build a binned representation of foils volume
    std::vector<double> binnedFoilsVolume();
    std::vector<double> weightedBinnedFoilsVolume();

    // methods to extract foil, position and time w.r.t. the chosen algorithm
    int getFlatRndFoil() ;
    int getVolumeRndFoil() ;
    int getExpoRndFoil() ;
    int getVolumeAndExpoRndFoil() ;
    CLHEP::Hep3Vector getFlatRndPos(TargetFoil const& theFoil) ;
    double getFlatRndTime() ;
    double getLimitedExpRndTime() ;
    
  };
} // end namespace mu2e,

#endif

