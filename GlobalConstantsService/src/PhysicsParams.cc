//
// Some physical parameters
//
//

// C++ includes
#include <algorithm>

// CLHEP includes
#include "CLHEP/Units/PhysicalConstants.h"

// Mu2e include files
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "Offline/ConfigTools/inc/SimpleConfig.hh"

// Framework includes
#include "cetlib/pow.h"
#include "cetlib_except/exception.h"

namespace mu2e {

  PhysicsParams::PhysicsParams( SimpleConfig const& config ) :
    _protonEnergy(0.),
    _protonKE(0.),
    _protonMomentum(0.)
  {

    const int verbosity = config.getInt("physicsParams.verbosityLevel",0);

    // Load physics constants
    const double alpha    = 1/config.getDouble("physicsParams.invFineStructureConstant");
    const double muonMass = config.getDouble("physicsParams.muonMass");
    const double pMass    = config.getDouble("physicsParams.protonMass");

    // Proton information
    _protonKE       = config.getDouble("physicsParams.protonKE");
    _protonEnergy   = _protonKE + pMass;
    _protonMomentum = std::sqrt( cet::diff_of_squares( _protonEnergy, pMass ) );

    _nominalDRPeriod = config.getDouble("physicsParams.nominalDRPeriod");
    _nominalDAQTicks  = config.getInt("physicsParams.nominalDAQTicks");
    _nominalDAQOnSpillBins  = config.getInt("physicsParams.nominalDAQOnSpillBins");
    _nominalDAQFrequency = config.getDouble("physicsParams.nominalDAQFrequency");

    std::vector<int> tmpPDGId;
    config.getVectorInt("physicsParams.freeLifePDGId", tmpPDGId, tmpPDGId);
    std::vector<double> tmpLifetime;
    config.getVectorDouble("physicsParams.freeLifetime", tmpLifetime, tmpLifetime);
    if(tmpPDGId.size() != tmpLifetime.size()) {
      throw cet::exception("BADINPUTS")
        <<"physicsParams.freeLifePDGId and physicsParams.freeLifetime must have the same number of entries\n";
    }
    for(unsigned i=0; i<tmpPDGId.size(); ++i) {
      freeLifetime_[PDGCode::type(tmpPDGId[i])] = tmpLifetime[i];
      // assume CPT
      if(tmpPDGId[i] < 0) {
        std::cerr<<"PhysicsParams() WARNING: ignoring PDG code sign for physicsParams.freeLifePDGId entry "
                 <<tmpPDGId[i]<<std::endl;
      }
      freeLifetime_[PDGCode::type(-tmpPDGId[i])] = tmpLifetime[i];
    }

    // Specify chosen stopping target; default is "Al"
    _chosenStoppingTargetMaterial = config.getString("physicsParams.chosenStoppingTargetMaterial","Al");

    // Load available stopping targets
    config.getVectorString("physicsParams.stoppingTargets",_allowedTargetMaterials);

    // Check if chosen stopping target is available; throw if not
    checkMaterial( _chosenStoppingTargetMaterial );

    // Loop over available data to fill maps
    std::string name;
    for ( const auto& material : _allowedTargetMaterials ) {

      // Load relevant atomic constants
      _atomicMass[material]     = config.getDouble("physicsParams."+material+".atomicMass")*CLHEP::amu_c2;
      _atomicNumber[material]   = config.getInt   ("physicsParams."+material+".atomicNumber");
      _decayTime[material]      = config.getDouble("physicsParams."+material+".decayTime" );
      _decayFraction[material]  = config.getDouble("physicsParams."+material+".decayFraction" );
      _atomicRadius[material]	= config.getDouble("physicsParams."+material+".atomicRadius");

      // Calculate approx. binding energy
      _approxBindingEnergy[material] = muonMass*cet::square(alpha*_atomicNumber[material])/2.;

      // Load muon energy; default to "mumass - Eb(approx.)" if not available
      const std::string muonEnergyString = "physicsParams."+material+".muonEnergy";
      bool goodMuonEnergy = config.hasName( muonEnergyString );
      if( goodMuonEnergy ) {
        _muonEnergy[material] = config.getDouble(muonEnergyString);
        _bindingEnergy[material] = muonMass - _muonEnergy[material];
      } else {
        _muonEnergy[material] = muonMass-_approxBindingEnergy[material] ;
      }

      // Calculate endpoint energy
      _endpointEnergy[material] = _muonEnergy[material] -cet::square(_muonEnergy[material])/(2*_atomicMass[material]);

      if ( verbosity > 0 ) {
        std::cout << " Stopping target parameters: " << material << std::endl
                  << "       - Muon Energy    : " << _muonEnergy[material] ;
        if ( ! goodMuonEnergy ) std::cout << " (using m_mu - m_mu*(Z*alpha)^2/2 approx.) " ;
        std::cout << std::endl
                  << "       - Endpoint Energy: " << _endpointEnergy[material] << std::endl;
      }

      // Load Czarnecki constants
      _czarneckiCoefficient[material] = config.getDouble("physicsParams."+material+".czarneckiCoefficient" );

      // the following quantities are optional for a material

      name = "physicsParams."+material+".czarneckiCoefficients";
      if(config.hasName(name)) config.getVectorDouble(name,_czarneckiCoefficients[material]);

      // Load capture product rates
      name = "physicsParams."+material+".capture.protonRate";
      if(config.hasName(name)) _captureProtonRate[material] = config.getDouble(name);
      name = "physicsParams."+material+".capture.deuteronRate";
      if(config.hasName(name)) _captureDeuteronRate[material] = config.getDouble(name);
      name = "physicsParams."+material+".capture.neutronRate";
      if(config.hasName(name)) _captureNeutronRate[material] = config.getDouble(name);
      name = "physicsParams."+material+".capture.photonRate";
      if(config.hasName(name)) _capturePhotonRate[material] = config.getDouble(name);

      // Load capture gamma rays
      name = "physicsParams."+material+".capture.photon.1809keV.energy";
      if(config.hasName(name)) _1809keVGammaEnergy[material] = config.getDouble(name);
      name = "physicsParams."+material+".capture.photon.1809keV.intensity";
      if(config.hasName(name)) _1809keVGammaIntensity[material] = config.getDouble(name);

      //Ce+ endpoint
      name = "physicsParams."+material+".ePlusEndpointEnergy";
      if(config.hasName(name)) _ePlusEndpointEnergy[material] = config.getDouble(name);

      name = "physicsParams."+material+".RMCbindingEnergyFit";
      if(config.hasName(name)) _RMCbindingEnergyFit[material] = config.getDouble(name);
      name = "physicsParams."+material+".RMCrecoilEnergyFit";
      if(config.hasName(name)) _RMCrecoilEnergyFit[material]  = config.getDouble(name);
      name = "physicsParams."+material+".RMCdeltaMassFit";
      if(config.hasName(name)) _RMCdeltaMassFit[material]     = config.getDouble(name);

    }

    // Load Shanker constants
    config.getVectorDouble("physicsParams.shankerDcoefficients", _shankerDcoefficients, this->getShankerNcoeffs() );
    config.getVectorDouble("physicsParams.shankerEcoefficients", _shankerEcoefficients, this->getShankerNcoeffs() );
    config.getVectorDouble("physicsParams.shankerFcoefficients", _shankerFcoefficients, this->getShankerNcoeffs() );

  }

  //================================================================
  double PhysicsParams::getParticleLifetime(PDGCode::type pdgId) const {
    FreeLifeMap::const_iterator i = freeLifetime_.find(pdgId);
    if(i == freeLifetime_.end()) {
      throw cet::exception("BADINPUTS")
        <<"PhysicsParams::getParticleLifetime(): no information for pdgId = "<<pdgId<<"\n";
    }
    return i->second;
  }

  //================================================================
  double PhysicsParams::doubleOrThrow( std::map<targetMat, double> const& tmap,
                                       targetMat const& material, std::string const& parameter) const {
      const std::string allowedMaterial = checkMaterial( material );
      auto it = tmap.find(allowedMaterial);
      if (it==tmap.end()) {
        throw cet::exception("PHYSICSPARAMS_MISSING_PAR")
          <<"PhysicsParams: request for missing parameter "<<parameter<<" for material "<<material<< "\n";
      }
      return it->second;
  }

}
