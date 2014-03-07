//
// StrawPhysics collects the electronics response behavior of a Mu2e straw in
// several functions.
//
// $Id: StrawPhysics.cc,v 1.6 2014/03/07 19:49:48 brownd Exp $
// $Author: brownd $
// $Date: 2014/03/07 19:49:48 $
//
// Original author David Brown, LBNL
//
#include "TrackerConditions/inc/StrawPhysics.hh"
#include "cetlib/exception.h"
#include <math.h>
#include <algorithm>

using namespace std;
namespace mu2e {
  StrawPhysics::StrawPhysics(fhicl::ParameterSet const& pset) :
    _EIonize(pset.get<double>("EnergyPerIonization",200.0e-6)), // 100% Ar is between 27 ev/ionization and 100 ev/ionization, not sure what model G4 uses, also should use Ar/CO2 FIXME!!
    _QIonize(pset.get<double>("ChargePerIonization",1.6e-7)), // e, pC
    _gasgain(pset.get<double>("GasGain",3.0e4)),
    _attlen{pset.get<double>("ShortAttentuationLength",50.0),pset.get<double>("LongAttentuationLength",27000.0)}, // from ATLAS TRT measurement
    _longfrac(pset.get<double>("LongAttenuationFraction",0.92)),
    _vdrift(pset.get<double>("DriftVelocity",0.05)), // mm/nsec
    _drifterr(pset.get<double>("DriftTimeError",1.5)), // nsec
    _vprop(pset.get<double>("PropagationVelocity",273.0)), //mm/nsec
    _vdisp(pset.get<double>("PropagationVelocityDispersion",0.01)) //1/nsec
  {}

  StrawPhysics::~StrawPhysics() {}
// for now, simple linear models of the effects.  In future, we want to
// include non-linear effects and angular dependence (ie Lorentz effects)
  double StrawPhysics::strawGain(double ddist, double phi) const { 
    return _gasgain;
  }

  double StrawPhysics::driftDistanceToTime(double ddist, double phi) const {
    return ddist/_vdrift;
  }

  double StrawPhysics::driftTimeSpread(double ddist, double phi) const {
    return _drifterr;
  }
  
  double StrawPhysics::propagationAttenuation(double wdist) const {
    return  (1.0-_longfrac)*exp(-wdist/_attlen[0]) +
      _longfrac*exp(-wdist/_attlen[1]);
  }

  double StrawPhysics::propagationTime(double wdist) const {
    return  wdist/_vprop;
  }

  double StrawPhysics::ionizationCharge(double ionizationEnergy) const {
    return _QIonize*ionizationEnergy/_EIonize; }

  double StrawPhysics::ionizationEnergy(double ionizationCharge) const {
    return _EIonize*ionizationCharge/_QIonize; }

  unsigned StrawPhysics::nIonization(double charge) const {
    return static_cast<unsigned>(std::max(1,static_cast<int>(ceil(charge/_QIonize))));
  }

}

