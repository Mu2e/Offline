#ifndef GlobalConstantsService_PhysicsParams_hh
#define GlobalConstantsService_PhysicsParams_hh
//
// Physical parameters used in Mu2e.
//
// NOTE: The muon and proton masses should come directly from the PDT;
// but I do not think those can be loaded when calculating the
// constants returned below as PhysicsParams and ParticleDataList are
// independent template arguments to GlobalConstantsService.  Because
// of this, the muon and proton masses must be specified in the
// SimpleConfig file which is passed as an argument to the
// constructor.  To discourage further use of these hard-coded mass
// values, no getProtonMass/getMuonMass accessors are provided.
// Thoses masses can be directly accessed later on through the
// GlobalConstantsHandle<ParticleDataList> construct.
//
// Author: Kyle Knoepfel
//

// C++ includes.
#include <map>
#include <vector>

// Mu2e includes.
#include "Offline/Mu2eInterfaces/inc/ConditionsEntity.hh"

// The use of this header  does not introduce a library dependency.
#include "Offline/DataProducts/inc/PDGCode.hh"

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

    double   getNominalDRPeriod() const { return _nominalDRPeriod; } // ns
    int      getNominalDAQTicks() const { return _nominalDAQTicks; } // counts, DAQ ticks per DRperiod
    int      getNominalDAQOnSpillBins() const { return _nominalDAQOnSpillBins; } // how often to 68->67
    double   getNominalDAQFrequency() const { return _nominalDAQFrequency; } // MHz
    double   getNominalDAQClockTick() const { return 1000.0/_nominalDAQFrequency; } // ns

    // Lifetimes of free (not stopped) particles.  We provide them
    // here because values coming from other sources are not accurate
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

    double getAtomicRadius  (targetMat material = "") const {
      const std::string allowedMaterial = checkMaterial( material );
      return _atomicRadius.find(allowedMaterial)->second;
    }

    double   getApproxEb      (targetMat material = "") const {
      const std::string allowedMaterial = checkMaterial( material );
      return _approxBindingEnergy.find(allowedMaterial)->second;
    }

    double   getEb            (targetMat material = "") const {
      return doubleOrThrow(_bindingEnergy,material,"BindingEnergy");
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
      auto it = _czarneckiCoefficients.find(allowedMaterial);
      if(it == _czarneckiCoefficients.end()) {
        throw cet::exception("PHYSICSPARAMS_MISSING_COEF")
          <<"PhysicsParams: request for missing Czarnecki coefficients for material "<<material<< "\n";
      }
      return it->second;
    }

    std::size_t getShankerNcoeffs() const { return _shankerNcoeffs; }
    const std::vector<double>& getShankerDcoefficients() const { return _shankerDcoefficients; }
    const std::vector<double>& getShankerEcoefficients() const { return _shankerEcoefficients; }
    const std::vector<double>& getShankerFcoefficients() const { return _shankerFcoefficients; }

    double   getCaptureProtonRate     (targetMat material = "") const {
      return doubleOrThrow(_captureProtonRate,material,"captureProtonRate");
    }
    double   getCaptureDeuteronRate     (targetMat material = "") const {
      return doubleOrThrow(_captureDeuteronRate,material,"captureDeuteronRate");
    }
    double   getCaptureNeutronRate     (targetMat material = "") const {
      return doubleOrThrow(_captureNeutronRate,material,"captureNeutronRate");
    }
    double   getCapturePhotonRate     (targetMat material = "") const {
      return doubleOrThrow(_capturePhotonRate,material,"capturePhotonRate");
    }

    double   get1809keVGammaEnergy     (targetMat material = "") const {
      return doubleOrThrow(_1809keVGammaEnergy,material,"1809keVGammaEnergy");
    }
    double   get1809keVGammaIntensity     (targetMat material = "") const {
      return doubleOrThrow(_1809keVGammaIntensity,material,"1809keVGammaIntensity");
    }

    double   getePlusEndpointEnergy(targetMat material = "") const {
      return doubleOrThrow(_ePlusEndpointEnergy,material,"ePlusEndpointEnergy");
    }

    double getRMCbindingEnergyFit(targetMat material = "")  const {
      return doubleOrThrow(_RMCbindingEnergyFit,material,"RMCbindingEnergyFit");
    }
    double getRMCrecoilEnergyFit(targetMat material = "") const {
      return doubleOrThrow(_RMCrecoilEnergyFit,material,"RMCrecoilEnergyFit");
    }
    double getRMCdeltaMassFit(targetMat material = "") const {
      return doubleOrThrow(_RMCdeltaMassFit,material,"RMCdeltaMassFit");
    }

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

    double _nominalDRPeriod;
    int    _nominalDAQTicks;
    int    _nominalDAQOnSpillBins;
    double _nominalDAQFrequency;

    typedef std::map<PDGCode::type, double> FreeLifeMap;
    FreeLifeMap freeLifetime_;

    std::map<targetMat,double>   _decayTime;
    std::map<targetMat,double>   _decayFraction;
    std::map<targetMat,double>   _atomicMass;
    std::map<targetMat,unsigned> _atomicNumber;
    std::map<targetMat,double>	 _atomicRadius;
    std::map<targetMat,double>   _approxBindingEnergy;
    std::map<targetMat,double>   _bindingEnergy;
    std::map<targetMat,double>   _muonEnergy;
    std::map<targetMat,double>   _endpointEnergy;
    std::map<targetMat,double>   _ePlusEndpointEnergy;

    std::map<targetMat,double> _czarneckiCoefficient;
    std::map<targetMat,std::vector<double>> _czarneckiCoefficients;

    const std::size_t _shankerNcoeffs = 4;
    std::vector<double> _shankerDcoefficients;
    std::vector<double> _shankerEcoefficients;
    std::vector<double> _shankerFcoefficients;

    std::map<targetMat, double> _captureProtonRate;
    std::map<targetMat, double> _captureDeuteronRate;
    std::map<targetMat, double> _captureNeutronRate;
    std::map<targetMat, double> _capturePhotonRate;

    std::map<targetMat, double> _1809keVGammaEnergy;
    std::map<targetMat, double> _1809keVGammaIntensity;

    std::map<targetMat, double> _RMCbindingEnergyFit;
    std::map<targetMat, double> _RMCrecoilEnergyFit;
    std::map<targetMat, double> _RMCdeltaMassFit;


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

    // check if parameter is available for this material, throw if not
    double doubleOrThrow(std::map<targetMat, double> const& tmap,
                         targetMat const& material, std::string const& parameter) const;

  };

}

#endif /* GlobalConstantsService_PhysicsParams_hh */
