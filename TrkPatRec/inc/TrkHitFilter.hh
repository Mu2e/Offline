//
// Simple class to filter hits based on residual,
//
//
// struct for outlier search tuple
#ifndef TrkHitFilter_hh
#define TrkHitFilter_hh
#include "CLHEP/Vector/ThreeVector.h"
#include "Rtypes.h"

namespace mu2e {
  struct TrkHitFilter {
    CLHEP::Hep3Vector _pos;
    Float_t _doca;
  };
}
#endif

