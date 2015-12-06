//
// Simple accessor to Kalman fit
//
// $Id: KalFitResult.hh,v 1.4 2012/12/05 21:14:14 tassiell Exp $
// $Author: tassiell $ 
// $Date: 2012/12/05 21:14:14 $
//
#ifndef KalFitResult_HH
#define KalFitResult_HH

#include "BTrk/BaBar/BaBar.hh"
// KalFit objects
#include "RecoDataProducts/inc/TrkFitDirection.hh"
#include "Mu2eBTrk/inc/TrkStrawHit.hh"
#include "TrkReco/inc/TrkDef.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/TrkBase/TrkParticle.hh"

// C++

namespace mu2e 
{
// struct defining the Kalman fit inputs and output
  struct KalFitResult {
    KalFitResult(): _tdef(0) ,_krep(0)   {}

// must initialize with a TrkDef
    KalFitResult(TrkDef const* tdef) : _tdef(tdef) ,_krep(0) {}
//-----------------------------------------------------------------------------
// KalFitResult doesn't own any pointers, '_krep' is handled in the pattern 
// recognition modules, so, in principle, no need to delete it here
// otherwise, a deletion of the list of KalFitResults (a data product) resutls in a crash
//-----------------------------------------------------------------------------
    ~KalFitResult() {}

    void    removeFailed() { if(_fit.failure())deleteTrack(); }
    void    fit() { _krep->fit(); }
    void    deleteTrack();
    KalRep* stealTrack() { KalRep* retval = _krep; _krep=0; _hits.clear(); return retval; }
// data payload
    TrkDef const*                _tdef; // original track definition on which this is based
    std::vector<DetIntersection> _detinter; // material intersections, used by the KalRep
  };
} 

#endif
