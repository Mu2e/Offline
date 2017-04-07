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
#include "BTrkData/inc/TrkStrawHit.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/TrkBase/TrkParticle.hh"
#include "RecoDataProducts/inc/Doublet.hh"

#include "CalPatRec/inc/TrkDefHack.hh"

#include <iostream>
// C++

namespace mu2e {
// struct defining the Kalman fit inputs and output
  struct KalFitResult {
//-----------------------------------------------------------------------------
// an internal CalPatRec data structure
// KalFitResult doesn't own any pointers, '_krep' is handled in the pattern 
// recognition modules, so, in principle, no need to delete it here
// otherwise, a deletion of the list of KalFitResults (a data product) resutls in a crash
//-----------------------------------------------------------------------------
    TrkDefHack const*                 _tdef;      // original track definition on which this is based
    const StrawHitCollection*         _shcol;     // 
    const StrawHitPositionCollection* _shpos;     //
    const StrawHitFlagCollection*     _shfcol;    //
    KalRep*                           _krep;      // Kalman rep, owned by the collection
    std::vector<TrkStrawHit*>         _hits;      // straw hits, owned by the KalRep, //
                                                  // now should not be needed as KalRep provides it
					          // however, it is still in use
    std::vector<DetIntersection> _detinter;       // material intersections, used by the KalRep
    TrkErrCode                   _fit;            // error code from last fit
    unsigned                     _nt0iter;        // number of times t0 was iterated
    unsigned                     _nweediter;      // number of iterations on hit weeding
    unsigned                     _nunweediter;    // number of iterations on hit unweeding
    unsigned                     _ninter;         // number of iterations on material intersections
    std::vector <Doublet>        _listOfDoublets; // list of hist multiplets
//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
    KalFitResult();
    KalFitResult(const TrkDefHack* hdef);
    ~KalFitResult();

    void    removeFailed() { if(_fit.failure()) deleteTrack(); }
    void    fit         () { if(_fit.success()) _fit = _krep->fit(); }
    void    deleteTrack ();
    
    KalRep* stealTrack             () { KalRep* retval = _krep; _krep=0; _hits.clear(); return retval; }
    const StrawHitCollection* shcol() { return _shcol; }
  };

  typedef std::vector<KalFitResult> KalFitResultCollection;
} 

#endif
