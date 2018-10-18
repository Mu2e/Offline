#ifndef GlobalConstantsService_PhysicsParams_hh
#define GlobalConstantsService_PhysicsParams_hh
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
#include <map>
#include <vector>

// Mu2e includes.
#include "Mu2eInterfaces/inc/ConditionsEntity.hh"

// The use of this header  does not introduce a library dependency.
#include "DataProducts/inc/PDGCode.hh"

// Framework includes
#include "cetlib_except/exception.h"

namespace mu2e
{
  class SimpleConfig;

  struct PhysicsParams: virtual public ConditionsEntity {

    typedef std::string targetMat;

    // Proton parameters
    double   getProtonEnergy  () const { return _protonEnergy;   }
    double   getProtonKE      () const { return _protonKE;       }
    double   getProtonMomentum() const { return _protonMomentum; }

    // Lifetimes of free (not stopped) particles.  We provide them
    // here because values coming from HepPDT are not accurate, and
    // are in wrong units.
    double getParticleLifetime(PDGCode::type pdgId) const;

    // Muon parameters
    double   getDecayTime     (targetMat material = "") const {
      const std::string allowedMaterial = checkMaterial( material );
      return _decayTime.find(allowedMaterial)->second;
    }

    double   getDecayFraction (targetMat material = "") const {
      const std::string allowedMaterial = checkMaterial( material );
      return _decayFraction.find(allowedMaterial)->second;
    }

    double   getAtomicMass    (targetMat material = "") const {
      const std::string allowedMaterial = checkMaterial( material );
      return _atomicMass.find(allowedMaterial)->second;
    }

    unsigned getAtomicNumber  (targetMat material = "") const {
      const std::string allowedMaterial = checkMaterial( material );
      return _atomicNumber.find(allowedMaterial)->second;
    }

    double   getApproxEb      (targetMat material = "") const {
      const std::string allowedMaterial = checkMaterial( material );
      return _approxBindingEnergy.find(allowedMaterial)->second;
    }

    double   getEb            (targetMat material = "") const {
      const std::string allowedMaterial = checkMaterial( material );
      return _bindingEnergy.find(allowedMaterial)->second;
    }

    double   getMuonEnergy    (targetMat material = "") const {
      const std::string allowedMaterial = checkMaterial( material );
      return _muonEnergy.find(allowedMaterial)->second;
    }

    double   getEndpointEnergy(targetMat material = "") const {
      const std::string allowedMaterial = checkMaterial( material );
      return _endpointEnergy.find(allowedMaterial)->second;
    }

    targetMat getStoppingTargetMaterial() const {
      return _chosenStoppingTargetMaterial;
    }

    // Return Czarnecki/Shanker coefficients
    double getCzarneckiCoefficient (targetMat material = "") const {
      const std::string allowedMaterial = checkMaterial( material );
      return _czarneckiCoefficient.find(allowedMaterial)->second;
    }

    const std::vector<double>& getCzarneckiCoefficients(targetMat material = "") const {
      const std::string allowedMaterial = checkMaterial( material );
      return _czarneckiCoefficients.find(allowedMaterial)->second;
    }

    std::size_t getShankerNcoeffs() const { return _shankerNcoeffs; }
    const std::vector<double>& getShankerDcoefficients() const { return _shankerDcoefficients; }
    const std::vector<double>& getShankerEcoefficients() const { return _shankerEcoefficients; }
    const std::vector<double>& getShankerFcoefficients() const { return _shankerFcoefficients; }

    PhysicsParams( SimpleConfig const& config );

    // Accept compiler generated:
    // copy c'tor, d'tor and assignment operator.

  private:

    // We want to discourage multi-phase construction.
    PhysicsParams ();

    targetMat _chosenStoppingTargetMaterial;
    std::vector<targetMat> _allowedTargetMaterials;

    double _protonEnergy;
    double _protonKE;
    double _protonMomentum;

    typedef std::map<PDGCode::type, double> FreeLifeMap;
    FreeLifeMap freeLifetime_;

    std::map<targetMat,double>   _decayTime;
    std::map<targetMat,double>   _decayFraction;
    std::map<targetMat,double>   _atomicMass;
    std::map<targetMat,unsigned> _atomicNumber;
    std::map<targetMat,double>   _approxBindingEnergy;
    std::map<targetMat,double>   _bindingEnergy;
    std::map<targetMat,double>   _muonEnergy;
    std::map<targetMat,double>   _endpointEnergy;

    std::map<targetMat,double> _czarneckiCoefficient;
    std::map<targetMat,std::vector<double>> _czarneckiCoefficients;

    const std::size_t _shankerNcoeffs = 4;
    std::vector<double> _shankerDcoefficients;
    std::vector<double> _shankerEcoefficients;
    std::vector<double> _shankerFcoefficients;

    inline targetMat checkMaterial( const targetMat& material ) const {
      if ( material.empty() ) return _chosenStoppingTargetMaterial;

      // Throw if not allowed
      else if ( std::find( _allowedTargetMaterials.begin(),
                           _allowedTargetMaterials.end(),
                           material ) == _allowedTargetMaterials.end() )
        throw cet::exception("StoppingTargetMaterial")
          << __func__ << " " << material << " not an allowed stopping target!\n" ;

      else return material;
    }

  };

}

#endif /* GlobalConstantsService_PhysicsParams_hh */
