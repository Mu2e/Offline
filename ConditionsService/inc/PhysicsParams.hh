#ifndef ConditionsService_PhysicsParams_hh
#define ConditionsService_PhysicsParams_hh
//
// Physical parameters used in Mu2e.  
//
// NOTE: The muon and proton masses should come directly from the PDT;
// but I do not think those can be loaded when calculating the
// constants returned below as PhysicsParams and ParticleDataTable are
// independent template arguments to GlobalConstantsService.  Because
// of this, the muon and proton masses must be specified in the
// SimpleConfig file which is passed as an argument to the
// constructor.  To discourage further use of these hard-coded mass
// values, no getProtonMass/getMuonMass accessors are provided.
// Thoses masses can be directly accessed later on through the
// GlobalConstantsHandle<ParticleDataTable> construct.
//
// Author: Kyle Knoepfel
//

// C++ includes.
#include <iostream>
#include <vector>

// Mu2e includes.
#include "Mu2eInterfaces/inc/ConditionsEntity.hh"

namespace mu2e
{
  class SimpleConfig;

  struct PhysicsParams: virtual public ConditionsEntity {

    // Proton parameters
    double   getProtonEnergy()   const { return _protonEnergy;   }
    double   getProtonKE()       const { return _protonKE;       }
    double   getProtonMomentum() const { return _protonMomentum; }
    
    // Muon parameters
    double   getDecayTime()      const { return _decayTime;      }
    double   getAtomicMass()     const { return _atomicMass;     }
    unsigned getAtomicNumber()   const { return _atomicNumber;   }

    double   getApproxEb()       const { return _approxBindingEnergy; }
    double   getEb()             const { return _bindingEnergy; }
    double   getMuonEnergy()     const { return _muonEnergy;     }
    double   getEndpointEnergy() const { return _endpointEnergy; }

    std::string getStoppingTarget() const { return _chosenStoppingTarget; }

    // Return Czarnecki/Shanker coefficients
    double getCzarneckiCoefficient() const { return _czarneckiCoefficient; }
    const std::vector<double>& getCzarneckiCoefficients() const { return _czarneckiCoefficients; }

    size_t getShankerNcoeffs() const { return _shankerNcoeffs; }
    const std::vector<double>& getShankerDcoefficients() const { return _shankerDcoefficients; }
    const std::vector<double>& getShankerEcoefficients() const { return _shankerEcoefficients; }
    const std::vector<double>& getShankerFcoefficients() const { return _shankerFcoefficients; }

    PhysicsParams( SimpleConfig const& config );

    // Accept compiler generated:
    // copy c'tor, d'tor and assignment operator.

  private:

    // We want to discourage multi-phase construction.
    PhysicsParams ();

    std::string _chosenStoppingTarget;
    std::vector<std::string> _allowedTargets;

    double _protonEnergy;
    double _protonKE;
    double _protonMomentum;

    double _decayTime;
    double   _atomicMass; 
    unsigned _atomicNumber;
    double   _approxBindingEnergy;
    double   _bindingEnergy;
    double   _muonEnergy;
    double   _endpointEnergy;

    double _czarneckiCoefficient;
    std::vector<double> _czarneckiCoefficients;

    const size_t _shankerNcoeffs = 4;
    std::vector<double> _shankerDcoefficients;
    std::vector<double> _shankerEcoefficients;
    std::vector<double> _shankerFcoefficients;

  };

  // Shift left (printing) operator.
  inline std::ostream& operator<<(std::ostream& ost,
                                  const PhysicsParams& lw ){
    ost << "( "
        << lw.getDecayTime() << ", "
        << " )";

    return ost;
  }

}

#endif /* ConditionsService_PhysicsParams_hh */
