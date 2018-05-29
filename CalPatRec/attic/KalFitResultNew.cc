///////////////////////////////////////////////////////////////////////////////
// Struct to hold BaBar Kalman fit
//
// $Id: KalFitResultNew.cc,v 1.1 2012/08/31 23:21:02 brownd Exp $
// $Author: brownd $ 
// $Date: 2012/08/31 23:21:02 $
//
// the following has to come before other BaBar includes
///////////////////////////////////////////////////////////////////////////////

#include "BTrk/BaBar/BbrCollectionUtils.hh"
#include "CalPatRec/inc/KalFitResultNew.hh"

namespace mu2e {

//-----------------------------------------------------------------------------
  KalFitResultNew::KalFitResultNew() {
    helixTraj   = NULL;
    tpart       = TrkParticle::e_minus;        // 11
    fdir        = TrkFitDirection::downstream; // = 0
    shcol       = 0;
    shpos       = 0;
    shfcol      = 0;
    krep        = 0;
    fit         = TrkErrCode::fail;
    nt0iter     = 0;
    nweediter   = 0;
    nunweediter = 0;

    hitIndices  = new vector<StrawHitIndex>;

    static unsigned icount(0);

    while (++icount < 10) {
      std::cout << "Hi Pasha, KalFitResultNew is DEPRECATED, please work on it." 
		<< " This message will repeat " << 10-icount << " more times."
		<< " In each job. Cheers, Dave." << std::endl;
    }
  }

//-----------------------------------------------------------------------------
  KalFitResultNew::~KalFitResultNew() {
    if (helixTraj) delete helixTraj;
    delete hitIndices;
  }

//-----------------------------------------------------------------------------
  void KalFitResultNew::deleteTrack() {
    if(krep != NULL) {
      delete krep;
      krep = NULL; 
    }
  } 

//-----------------------------------------------------------------------------
  void KalFitResultNew::init() {
    deleteTrack();
    hitIndices->clear();
    missingHits.clear();
    //    doca.clear();
  } 

//-----------------------------------------------------------------------------
  KalRep*  KalFitResultNew::stealTrack() { 
    KalRep* retval = krep; 
    krep           = 0; 
    return retval; 
  }

}

