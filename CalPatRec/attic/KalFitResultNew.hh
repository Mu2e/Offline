//
// Simple accessor to Kalman fit
//
// $Id: $
// $Author: $ 
// $Date:  $
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
  struct MissingHit_t {
    StrawHitIndex  index;
    double         doca;
    double         dr;
  };

//-----------------------------------------------------------------------------
// struct defining the Kalman fit inputs and output
// an internal CalPatRec data structure
// KalFitResultNew doesn't own any pointers, '_krep' is handled in the pattern 
// recognition modules, so, in principle, no need to delete it here
// otherwise, a deletion of the list of KalFitResultNews (a data product) resutls in a crash
//-----------------------------------------------------------------------------
  struct KalFitResultNew {
    const art::Event*                 event;
    KalRep*                           krep;           // Kalman rep, owned by the collection
    const StrawHitCollection*         shcol;          // 
    const StrawHitPositionCollection* shpos;          //
    const StrawHitFlagCollection*     shfcol;         //
    std::string                       shDigiLabel;    // 
    TrkParticle                       tpart;
    TrkFitDirection                   fdir;
    const CaloCluster*                caloCluster;    //

    const HelixSeed*                  helixSeed;      //
    KalSeed*                          kalSeed;        // 
    TrkT0                             t0;             // estimate of the track t0
    HelixTraj*                        helixTraj;      // initial parameterization of the track
    std::vector<StrawHitIndex>*       hitIndices;     // list of hit indices, updates during the fit
    std::vector<StrawHitIndex>*       savedHits;      // list of hit indices, updates during the fit

    TrkErrCode                        fit;            // error code from last fit
    unsigned                          nt0iter;        // number of times t0 was iterated
    unsigned                          nweediter;      // number of iterations on hit weeding
    unsigned                          nunweediter;    // number of iterations on hit unweeding
    std::vector <Doublet>             listOfDoublets; // list of hist multiplets
    int                               nrescued;       // N rescued hits
    // std::vector<StrawHitIndex>        missingHits;    // used by findMissingHits and addHits
    // std::vector<double>               doca;           // used by findMissingHits and addHits
    std::vector<MissingHit_t>         missingHits; 
    int                               fitType;        // 0:seed 1:final
//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
    KalFitResultNew();
    ~KalFitResultNew();

    void    removeFailed() { if(fit.failure()) deleteTrack(); }
    void    deleteTrack ();
    KalRep* stealTrack  ();
    void    init        ();

    int                               nHelixHits     () { return helixSeed->hits().size(); }
    std::vector<StrawHitIndex>*       strawHitIndices() { return hitIndices; }

    //    const HelixTraj*                  helixTraj      () { return _helixTraj; }
  };

} 

#endif
