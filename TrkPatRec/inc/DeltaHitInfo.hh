//
// $Id: DeltaHitInfo.hh,v 1.1 2012/09/24 18:39:55 brownd Exp $
// $Author: brownd $ 
// $Date: 2012/09/24 18:39:55 $
//
// struct for hit diagnostics
#ifndef DeltaHitInfo_hh
#define DeltaHitInfo_hh
//#include "BaBar/BaBar.hh"
#include "TrkPatRec/inc/StrawHitInfo.hh"
namespace mu2e {
typedef StrawHitInfo DeltaHitInfo;
/*  struct DeltaHitInfo : StrawHitInfo {
    Int_t _hflag;
    Float_t _hgd; // MVA output of generalized distance
    Float_t _dphi, _drho, _dt; // MVA inputs
    virtual ~DeltaHitInfo();
// root
    ClassDef(DeltaHitInfo,1)
  };
  */
}
#endif
