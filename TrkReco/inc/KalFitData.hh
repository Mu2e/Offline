//
// Simple accessor to Kalman fit
//
//
#ifndef KalFitData_HH
#define KalFitData_HH

#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/TrkFitDirection.hh"
#include "RecoDataProducts/inc/HelixSeed.hh"
#include "RecoDataProducts/inc/KalSeed.hh"
#include "BTrkData/inc/Doublet.hh"

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
// KalFitData doesn't own any pointers, '_krep' is handled in the pattern 
// recognition modules, so, no need to delete it here
// otherwise, deletion of the list of KalFitData's (a data product) results in a crash
//-----------------------------------------------------------------------------
  struct KalFitData {
    struct Diag_t {
      unsigned  diskId;	
      unsigned  added;	
      double    depth;	
      double    dt;	
      double    trkPath;
      double    energy;	
      double    doca;   
    };
    
    const art::Event*                 event;
    KalRep*                           krep;           // Kalman rep, owned by the collection
    const ComboHitCollection*         chcol;          // 

    const StrawHitFlagCollection*     shfcol;         //
    std::string                       shDigiLabel;    // 

    TrkFitDirection                   fdir;
    const CaloCluster*                caloCluster;    //
    const CaloClusterCollection*      caloClusterCol;    //

    const HelixSeed*                  helixSeed;      //
    const KalSeed*                    kalSeed;        // 
    HelixTraj*                        helixTraj;      // initial parameterization of the track
    unsigned                          nweediter;      // number of iterations on hit weeding
    unsigned                          nweedtchiter;   // number of iterations on TrkCaloHit weeding
    std::vector<MissingHit_t>         missingHits; 
    int                               fitType;        // 0:seed 1:final
    
    Diag_t                            diag;
//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
    KalFitData();
    ~KalFitData();

    void    deleteTrack ();
    KalRep* stealTrack  ();
    void    init        ();
  };

} 

#endif
