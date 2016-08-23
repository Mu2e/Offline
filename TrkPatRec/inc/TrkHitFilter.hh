//
// Simple class to filter hits based on residual, 
//
// $Id: TrkHitFilter.hh,v 1.3 2014/05/05 22:25:56 brownd Exp $
// $Author: brownd $ 
// $Date: 2014/05/05 22:25:56 $
//
// struct for outlier search tuple
#ifndef TrkHitFilter_hh
#define TrkHitFilter_hh
#include "BTrk/BaBar/BaBar.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include "Rtypes.h"

namespace mu2e {
  struct TrkHitFilter {
    CLHEP::Hep3Vector _pos;
    Float_t _doca;
  };
}
#endif

