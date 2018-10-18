#ifndef ConditionsService_TrackerCalibrationStructs_hh
#define ConditionsService_TrackerCalibrationStructs_hh
//
// Parameters for tracker calibrations.
//
// Original author David Brown
//

// C++ includes.
#include <iostream>

// Mu2e includes.
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e
{
  class SimpleConfig;
  class StrawHit;
  class Straw;
// simple struct to hold output of timeToDistance function
  struct T2D {
    double _rdrift;
    double _rdrifterr;
    double _vdrift; // local drift velocity at this radius, ie dr/dt(rdrift).
    T2D() : _rdrift(0.0), _rdrifterr(1.0), _vdrift(0.0) {}
  };

// simple struct to hold output of distanceToTime function
  struct D2T {
    double _tdrift;
    double _tdrifterr;
    double _vdrift; // local drift velocity at this radius, ie dr/dt(rdrift).
    D2T() : _tdrift(0.0),_tdrifterr(1.0), _vdrift(0.0) {}
  };

// simple struct to hold output of energyToAmplitude function
  struct E2A {
    double _ampl;
    double _amplerr;
    E2A() : _ampl(0.0), _amplerr(1.0) {}
  };

// simple struct to hold output of amplitudeToEnergy function
  struct A2E {
    double _edep;
    double _edeperr;
    A2E() : _edep(0.0),_edeperr(1.0) {}
  };

}

#endif
