//
// Simple accessor to Kalman fit
//
// $Id: KalFitResultNew.hh,v 1.4 2012/12/05 21:14:14 tassiell Exp $
// $Author: tassiell $ 
// $Date: 2012/12/05 21:14:14 $
//
#ifndef KalFitResultNew_HH
#define KalFitResultNew_HH

#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/TrkFitDirection.hh"
#include "RecoDataProducts/inc/HelixSeed.hh"
#include "RecoDataProducts/inc/KalSeed.hh"
#include "RecoDataProducts/inc/Doublet.hh"

#include "BTrk/TrkBase/TrkT0.hh"
#include "BTrk/BaBar/BaBar.hh"

#include "BTrkData/inc/TrkStrawHit.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/TrkBase/TrkParticle.hh"
#include "BTrk/TrkBase/HelixTraj.hh"

namespace art {
  class Event;
}

namespace mu2e {
//-----------------------------------------------------------------------------
// struct defining the Kalman fit inputs and output
// an internal CalPatRec data structure
// KalFitResultNew doesn't own any pointers, '_krep' is handled in the pattern 
// recognition modules, so, in principle, no need to delete it here
// otherwise, a deletion of the list of KalFitResultNews (a data product) resutls in a crash
//-----------------------------------------------------------------------------
  struct KalFitResultNew {
    const art::Event*                 _event;
    KalRep*                           _krep;           // Kalman rep, owned by the collection
    const StrawHitCollection*         _shcol;          // 
    const StrawHitPositionCollection* _shpos;          //
    const StrawHitFlagCollection*     _shfcol;         //
    std::string                       _shDigiLabel;    // 
    TrkParticle                       _tpart;
    TrkFitDirection                   _fdir;
    const CaloCluster*                _caloCluster;    //

    const HelixSeed*                  _helixSeed;      //
    KalSeed*                          _kalSeed;        // 
    TrkT0                             _t0;             // estimate of the track t0
    HelixTraj*                        _helixTraj;      // initial parameterization of the track
    std::vector<StrawHitIndex>*       _hitIndices;     // list of hit indices, updates during the fit
    std::vector<StrawHitIndex>*       _savedHits;      // list of hit indices, updates during the fit

    TrkErrCode                        _fit;            // error code from last fit
    unsigned                          _nt0iter;        // number of times t0 was iterated
    unsigned                          _nweediter;      // number of iterations on hit weeding
    unsigned                          _nunweediter;    // number of iterations on hit unweeding
    std::vector <Doublet>             _listOfDoublets; // list of hist multiplets
    int                               _nrescued;       // N rescued hits
    std::vector<StrawHitIndex>        _missingHits;    // used by findMissingHits and addHits
    std::vector<double>               _doca;           // used by findMissingHits and addHits
    int                               _fitType;        // 0:seed 1:final
//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
    KalFitResultNew();
    ~KalFitResultNew();

    void    removeFailed() { if(_fit.failure()) deleteTrack(); }
    void    deleteTrack ();
    KalRep* stealTrack  ();
    void    init        ();

    const StrawHitCollection*         shcol () { return _shcol;  }
    const StrawHitPositionCollection* shpos () { return _shpos ; }
    const StrawHitFlagCollection*     shfcol() { return _shfcol; }

    const CaloCluster*                caloCluster() { return _caloCluster; }

    int                               nHelixHits     () { return _helixSeed->hits().size(); }
    std::vector<StrawHitIndex>*       strawHitIndices() { return _hitIndices; }
    TrkT0&                            t0()              { return _t0; }

    const HelixTraj*                  helixTraj      () { return _helixTraj; }
  };

} 

#endif
