//
//
// Simple class to filter hits based on residual, 
//
// $Id: TrkHitFilter.hh,v 1.1 2011/09/27 21:58:01 mu2ecvs Exp $
// $Author: mu2ecvs $ 
// $Date: 2011/09/27 21:58:01 $
//
// struct for outlier search tuple
#include "BaBar/BaBar.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include "Rtypes.h"

namespace mu2e {
  struct TrkHitInfo {
    CLHEP::Hep3Vector _pos;
    Float_t _resid;
    Int_t _mcpdg;
    Int_t _mcgen;
    Int_t _mcproc;    
// root 
    ClassDef(TrkHitInfo,1)
  };
}
