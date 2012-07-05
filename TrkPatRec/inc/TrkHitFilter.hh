//
//
// Simple class to filter hits based on residual, 
//
// $Id: TrkHitFilter.hh,v 1.2 2012/07/05 21:38:53 brownd Exp $
// $Author: brownd $ 
// $Date: 2012/07/05 21:38:53 $
//
// struct for outlier search tuple
#ifndef TrkHitFilter_hh
#define TrkHitFilter_hh
#include "BaBar/BaBar.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include "Rtypes.h"

namespace mu2e {
  struct TrkHitFilter {
    CLHEP::Hep3Vector _pos;
    Float_t _doca;
    Int_t _mcpdg;
    Int_t _mcgen;
    Int_t _mcproc;    
// root 
    ClassDef(TrkHitFilter,1)
  };
}
#endif

