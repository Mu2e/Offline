#ifndef FOILPARTICLEGENERATOR_HH
#define FOILPARTICLEGENERATOR_HH

//
// Generate a generic particle coming from target foils.
// Position, time and foil number of generated particle 
// is extracted random from appropiate distributions
//
//
// For now this is limited to:
//  - Uniform over the targets.
//  - Uniform in time during the requested time interval.
//  - All of the above need to be improved at a later date.
//  - Limits on cos(theta) and phi but uniform within the range.
 

#include <memory>

// Framework includes
#include "FWCore/Services/interface/RandomNumberGeneratorService.h"

// Mu2e includes
#include "Mu2eUtilities/inc/PDGCode.hh"
#include "ToyDP/inc/ToyGenParticleCollection.hh"
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
    
    FoilParticleGenerator( PDGCode::type pdgId, 
                           GenId::enum_type genId);
    
    ~FoilParticleGenerator();

    void generateFromFoil(ToyGenParticleCollection& genParts, long nParticles = 1);

    void setRandomEngine( edm::RandomNumberGeneratorService::base_engine_t& engine );

  private:


    FoilParticleGenerator();

    PDGCode::type _pdgId;
    GenId::enum_type _genId;

    // Random numbers generators
    std::auto_ptr<CLHEP::RandFlat>     _randFlat;
    std::auto_ptr<CLHEP::RandExponential>     _randTime;
    std::auto_ptr<RandomUnitSphere>    _randomUnitSphere;
    std::auto_ptr<CLHEP::RandGeneral>  _shape;
    std::auto_ptr<CLHEP::RandGeneral>  _randFoils;

    //number of foils
    int _nfoils;

    //
    double _conversionEnergyAluminum;

    //Limits on energy range and number of bins representation of energy spectrum
    double _elow, _ehi;
    int _nbins;

    //Spatial and temporal variables
    double _czmin, _czmax, _phimin, _phimax, _tmin, _tmax, _p, _mass;

    //Range for the above
    double _dcz, _dphi, _dt;

    // Compute the value of the energy spectrum at given energy.
    double energySpectrum(double e);

    // Build a binned representation of the energy spectrum.
    std::vector<double> binnedEnergySpectrum();

    //Build a binned representation of foils volume
    std::vector<double> binnedFoilsVolume();


  };
} // end namespace mu2e,

#endif

