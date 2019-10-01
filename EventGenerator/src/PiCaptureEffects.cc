//
//
// $Id: PiCaptureEffects.cc,v 1.2 2012/07/26 19:01:01 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/07/26 19:01:01 $
//
// Original author: Gianni Onorato
//

// C++ includes
#include <iostream>

// Framework includes
#include "cetlib_except/exception.h"

// Mu2e includes
#include "EventGenerator/inc/PiCaptureEffects.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GeneralUtilities/inc/safeSqrt.hh"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

using namespace std;

namespace mu2e {

  PiCaptureEffects::PiCaptureEffects(CLHEP::HepRandomEngine& engine,
                                     double const probInternalConversion,
                                     double const elow,
                                     double const ehi,
                                     int const nbins):
    _probInternalConversion(probInternalConversion),
    _elow(elow),
    _ehi(ehi),
    _nbins (nbins),
    _randomUnitSphere{engine},
    _randFlat{engine},
    _spectrum{engine, &(binnedEnergySpectrum()[0]), _nbins},
    _internalFractionalSpectrum{engine, &(internalFractionalBinnedSpectrum()[0]), _nbins}
  {
    GlobalConstantsHandle<ParticleDataTable> pdt;
    // pick up particle mass
    const HepPDT::ParticleData& e_data = pdt->particle(PDGCode::e_minus).ref();
    _electMass = e_data.mass().value();
  }

  void PiCaptureEffects::defineOutput() {

    _e = _elow + _spectrum.fire() * (_ehi - _elow);
    _fract = _internalFractionalSpectrum.fire();
    _doPhoton = _randFlat.fire() > _probInternalConversion;
    _eElec = _e * _fract;
    _ePosit = _e * ( 1 - _fract);
    _elecMom = safeSqrt((_eElec*_eElec) - (_electMass*_electMass));
    _positMom = safeSqrt((_ePosit*_ePosit) - (_electMass*_electMass));

  }


  double PiCaptureEffects::internalFraction() const {
    return _fract;
  }

  bool PiCaptureEffects::doPhoton() const {
    return _doPhoton;
  }

  bool PiCaptureEffects::doElectron() const {
    return (_elecMom > 0 && !_doPhoton);
  }

  bool PiCaptureEffects::doPositron() const {
    return (_positMom > 0 && !_doPhoton);
  }


  GenParticle PiCaptureEffects::outputGamma(CLHEP::Hep3Vector pos, double time) {
    if (!_doPhoton) {
      throw cet::exception("MODEL") << "Internal conversion: you have to ask for pair production!\n";
    }
    CLHEP::HepLorentzVector mom (_randomUnitSphere.fire(_e), _e);

    GenParticle outGen(PDGCode::gamma, GenId::ExternalRPC, pos, mom, time);

    return outGen;
  }


  GenParticle PiCaptureEffects::outputElec(CLHEP::Hep3Vector pos, double time) {

    if (!doElectron()) {
      throw cet::exception("MODEL") << "Electron not created. Check with doElectron() before creating it.\n";
    }

    CLHEP::HepLorentzVector mom (_randomUnitSphere.fire(_elecMom), _elecMom);

    GenParticle outGen(PDGCode::e_minus, GenId::InternalRPC, pos, mom, time);

    return outGen;
  }

  GenParticle PiCaptureEffects::outputPosit(CLHEP::Hep3Vector pos, double time) {

    if (!doPositron()) {
      throw cet::exception("MODEL") << "Positron not created. Check with doElectron() before creating it.\n";
    }

    CLHEP::HepLorentzVector mom (_randomUnitSphere.fire(_positMom), _positMom);

    GenParticle outGen(PDGCode::e_plus, GenId::InternalRPC, pos, mom, time);

    return outGen;
  }


  // Photon energy spectrum as a continuous function.
  double PiCaptureEffects::energySpectrum(const double x) const
  {
    // Parameters from doc 665-v2
    constexpr double emax  = 138.2;
    constexpr double alpha =   2.691;
    constexpr double gamma =   1.051;
    constexpr double tau   =   8.043;
    constexpr double c0    =   2.741;
    constexpr double c1    =  -0.005;

    return pow(emax-x,alpha) * exp(-(emax-gamma*x)/tau) * (c0 + c1*x);

  } // PiCaptureEffects::energySpectrum

  // Compute a binned representation of the photon energy spectrum.
  std::vector<double> PiCaptureEffects::binnedEnergySpectrum() const{

    // Sanity check.
    if (_nbins <= 0) {
      throw cet::exception("RANGE")
        << "Nonsense PiCaptureEffectsGun.nbins requested="
        << _nbins
        << "\n";
    }

    // Bin width.
    double dE = (_ehi - _elow) / _nbins;

    // Vector to hold the binned representation of the energy spectrum.
    std::vector<double> spectrum;
    spectrum.reserve(_nbins);

    for (int ib=0; ib<_nbins; ib++) {
      double x = _elow+(ib+0.5) * dE;
      spectrum.push_back(energySpectrum(x));
    }

    return spectrum;
  } // PiCaptureEffects::binnedEnergySpectrum



  // Photon energy spectrum as a continuous function.
  double PiCaptureEffects::internalFractionalSpectrum(const double x) const
  {
    // the usual distribution for positron/electron fraction, normalized to unity
    return (9./7.)*(1. - (4/3)*x*(1 - x));
  } // PiCaptureEffects::internalFractionalSpectrum

  // Compute a binned representation of the photon energy spectrum.
  // This runs from 0 to 1; we'll take the photon energy and multiply.
  std::vector<double> PiCaptureEffects::internalFractionalBinnedSpectrum() const {

    // Sanity check.
    if (_nbins <= 0) {
      throw cet::exception("RANGE")
        << "Nonsense InternalPiCaptureEffectsGun.nbins requested="
        << _nbins
        << "\n";
    }

    // Bin width.
    double dy = 1./_nbins;

    // Vector to hold the binned representation of the energy spectrum.
    std::vector<double> spectrum;
    spectrum.reserve(_nbins);

    for (int ib=0; ib<_nbins; ib++) {
      double x = (ib+0.5) * dy;
      spectrum.push_back(internalFractionalSpectrum(x));
    }

    return spectrum;
  } // PiCaptureEffects::internalFractionalBinnedSpectrum

 }
