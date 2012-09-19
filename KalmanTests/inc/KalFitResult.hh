//
// Simple accessor to Kalman fit
//
// $Id: KalFitResult.hh,v 1.2 2012/09/19 20:17:37 brownd Exp $
// $Author: brownd $ 
// $Date: 2012/09/19 20:17:37 $
//
#ifndef KalFitResult_HH
#define KalFitResult_HH

#include "BaBar/BaBar.hh"
// KalFit objects
#include "KalmanTests/inc/TrkFitDirection.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"
#include "KalmanTests/inc/TrkDef.hh"
#include "KalmanTrack/KalRep.hh"
#include "TrkBase/TrkParticle.hh"
// C++

namespace mu2e 
{
// struct defining the Kalman fit inputs and output
  struct KalFitResult {
// must initialize with a TrkDef
    KalFitResult(TrkDef const& tdef) : _tdef(tdef) ,_krep(0), _fit(TrkErrCode::fail), _nt0iter(0), _nweediter(0), _nunweediter(0) {}
    ~KalFitResult() { delete _krep;}
    void removeFailed() { if(_fit.failure())deleteTrack(); }
    void fit() { if(_fit.success()) _fit = _krep->fit(); }
    void deleteTrack();
    KalRep* stealTrack() { KalRep* retval = _krep; _krep=0; _hits.clear(); return retval; }
// data payload
    TrkDef const& _tdef; // original track definition on which this is based
    KalRep* _krep; // Kalman rep, owned by the collection
    std::vector<TrkStrawHit*> _hits; // straw hits, owned by the KalRep
    std::vector<DetIntersection> _detinter; // material intersections, used by the KalRep
    TrkErrCode _fit; // error code from last fit
    unsigned _nt0iter; // number of times t0 was iterated
    unsigned _nweediter; // number of iterations on hit weeding
    unsigned _nunweediter; // number of iterations on hit unweeding
  };
} 

#endif
