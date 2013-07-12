//
// Some physical parameters
//
//

// CLHEP includes
#include "CLHEP/Units/PhysicalConstants.h"

// Mu2e include files
#include "ConditionsService/inc/PhysicsParams.hh"
#include "ConfigTools/inc/SimpleConfig.hh"

// Framework includes
#include "cetlib/exception.h"
#include "cetlib/pow.h"

namespace mu2e {

  PhysicsParams::PhysicsParams( SimpleConfig const& config ) :
    _decayTime(0.), 
    _atomicMass(0.), 
    _atomicNumber(0), 
    _approxBindingEnergy(0.), 
    _bindingEnergy(0.),
    _muonEnergy(0.),
    _endpointEnergy(0.) 
  {

    const int verbosity = config.getInt("physicsParams.verbosityLevel",0);
    
    // Load available stopping targets
    config.getVectorString("physicsParams.stoppingTargets",_allowedTargets);

    // Specify chosen stopping target; default is "Al"
    _chosenStoppingTarget = config.getString("physicsParams.chosenStoppingTarget","Al");

    // Check if chosen stopping target is available; throw if not
    if ( std::find( _allowedTargets.begin(),
                    _allowedTargets.end(),
                    _chosenStoppingTarget ) == _allowedTargets.end() )
      throw cet::exception("StoppingTargetMaterial")
        << __func__ << " " << _chosenStoppingTarget << " not an allowed stopping target!\n" ;

    // Load physics constants
    const double alpha    = 1/config.getDouble("physicsParams.invFineStructureConstant");
    const double muonMass = config.getDouble("physicsParams.muonMass");

    // Load relevant atomic constants
    _atomicMass     = config.getDouble("physicsParams."+_chosenStoppingTarget+".atomicMass")*CLHEP::amu_c2;
    _atomicNumber   = config.getInt   ("physicsParams."+_chosenStoppingTarget+".atomicNumber");
    _decayTime      = config.getDouble("physicsParams."+_chosenStoppingTarget+".decayTime" );

    // Calculate approx. binding energy
    _approxBindingEnergy = muonMass*cet::square(alpha*_atomicNumber)/2.;

    // Load muon energy; default to "mumass - Eb(approx.)" if not available
    const std::string muonEnergyString = "physicsParams."+_chosenStoppingTarget+".muonEnergy";
    const bool approxRequired = !config.hasName( muonEnergyString );
    _muonEnergy     = config.getDouble( muonEnergyString , muonMass-_approxBindingEnergy );
    
    // Calculate actual binding energy if available
    if ( !approxRequired ) _bindingEnergy = muonMass - _muonEnergy;

    // Calculate endpoint energy
    _endpointEnergy = _muonEnergy -cet::square(_muonEnergy)/(2*_atomicMass);

    if ( verbosity > 0 ) {
      std::cout << " Chosen stopping target : " << _chosenStoppingTarget << std::endl
                << "       - Muon Energy    : " << _muonEnergy ;
      if ( approxRequired ) std::cout << " (using m_mu - m_mu*(Z*alpha)^2/2 approx.) " ;
      std::cout << std::endl
                << "       - Endpoint Energy: " << _endpointEnergy << std::endl;
    }

    // Load Czarnecki constants
    _czarneckiCoefficient = config.getDouble("physicsParams."+_chosenStoppingTarget+".czarneckiCoefficient" );
    config.getVectorDouble("physicsParams."+_chosenStoppingTarget+".czarneckiCoefficients",_czarneckiCoefficients,std::vector<double>());

    // Load Shanker constants
    config.getVectorDouble("physicsParams.shankerDcoefficients", _shankerDcoefficients, this->getShankerNcoeffs() );
    config.getVectorDouble("physicsParams.shankerEcoefficients", _shankerEcoefficients, this->getShankerNcoeffs() );
    config.getVectorDouble("physicsParams.shankerFcoefficients", _shankerFcoefficients, this->getShankerNcoeffs() );

  }

}
