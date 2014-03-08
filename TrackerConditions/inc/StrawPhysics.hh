#ifndef TrackerConditions_StrawPhysics_hh
#define TrackerConditions_StrawPhysics_hh
//
// StrawPhysics collects the electronics response behavior of a Mu2e straw in
// several functions and parameters
//
// $Id: StrawPhysics.hh,v 1.2 2014/03/08 00:55:21 brownd Exp $
// $Author: brownd $
// $Date: 2014/03/08 00:55:21 $
//
// Original author David Brown, LBNL
//

// C++ includes
#include <iostream>
#include <array>
#include <vector>
#include <utility>

// Mu2e includes
#include "Mu2eInterfaces/inc/ConditionsEntity.hh"
#include "fhiclcpp/ParameterSet.h"
#include <array>
#include <vector>

namespace mu2e {
  class StrawPhysics : virtual public ConditionsEntity {
    public:
      // construct from parameters
      StrawPhysics(fhicl::ParameterSet const& pset);
      virtual ~StrawPhysics();
    // models.  Note these are different from the corresponding
    // functions used in reconstruction, as those can be wire-
    // dependent and have emergent properties
      double ionizationCharge(double ionizationEnergy) const;
      double ionizationEnergy(double ionizationCharge) const;
      unsigned nIonization(double charge) const; 
      double strawGain(double ddist, double phi) const;
      double driftDistanceToTime(double ddist, double phi) const;  // single cluster!
      double driftTimeSpread(double ddist, double phi) const; // single cluster!
      double propagationAttenuation(double wdist) const; 
      double propagationTime(double wdist) const;
      double velocityDispersion() const { return _vdisp; } 

    private:
      double _EIonize; // energy of each ionization (MeV)
      double _QIonize; // charge of a single ionization (=e, pC)
      double _gasgain; // avalanche gain
      // attenuation length of charge down the wire; note
      // there is a short and a long component, each with it's own amplitude
      double _attlen[2];
      double _longfrac;
      double _vprop; // (average) propagation velocity
      double _vdisp; // dispersion of propagation velocity (dv/dl)
    // parameters describing cluster DtoT
      std::vector<double> _cdpoly;
      std::vector<double> _cdsigmapoly;

  };
}
#endif

