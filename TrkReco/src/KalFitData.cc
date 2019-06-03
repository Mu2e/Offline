///////////////////////////////////////////////////////////////////////////////
// Struct to hold BaBar Kalman fit
//
// $Id: $
// $Author: $ 
// $Date: $
//
// the following has to come before other BaBar includes
///////////////////////////////////////////////////////////////////////////////

#include "BTrk/BaBar/BbrCollectionUtils.hh"
#include "TrkReco/inc/KalFitData.hh"

namespace mu2e {

//-----------------------------------------------------------------------------
  KalFitData::KalFitData() {
    helixTraj   = NULL;
    //    tpart       = TrkParticle::e_minus;        // 11
    fdir        = TrkFitDirection::downstream; // = 0
    chcol       = 0;
    // shpos       = 0;
    shfcol      = 0;
    krep        = 0;
    kalSeed     = 0;
    helixSeed   = 0;
    //    fit         = TrkErrCode::fail;
    //    nt0iter     = 0;
    nweediter   = 0;
    nweedtchiter   = 0;
    //    nunweediter = 0;

    // hitIndices  = new vector<StrawHitIndex>;
  }

//-----------------------------------------------------------------------------
  KalFitData::~KalFitData() {
    // if (helixTraj)   { delete helixTraj; helixTraj = 0;}
    // if (caloCluster) { delete caloCluster; caloCluster = 0;}
    // if (helixSeed)   { delete helixSeed; helixSeed  = 0;}
    // if (kalSeed)     { delete kalSeed; kalSeed = 0;}
    // if (shcol)       { delete shcol; shcol = 0;}
    // if (shfcol)      { delete shfcol; shfcol = 0;}
    // if (event)       { delete event; event = 0;}
    // if (krep)        { delete krep; krep = 0;}

    // delete hitIndices;
  }

//-----------------------------------------------------------------------------
  void KalFitData::deleteTrack() {
    if(krep != NULL) {
      delete krep;
      krep = NULL; 
    }
  } 

//-----------------------------------------------------------------------------
  void KalFitData::init() {
    deleteTrack();
    missingHits.clear();
    nweediter    = 0;
    nweedtchiter = 0;
    helixTraj    = NULL;
    
    diag.diskId   = 0;
    diag.added	  = 0;
    diag.depth    = -9999;
    diag.dt       = -9999;
    diag.trkPath  = -9999;
    diag.energy   = -9999;
    diag.doca     = -9999;   
}

//-----------------------------------------------------------------------------
  KalRep*  KalFitData::stealTrack() { 
    KalRep* retval = krep; 
    krep           = 0; 
    return retval; 
  }

}

