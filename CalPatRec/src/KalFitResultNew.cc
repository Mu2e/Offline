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
    _helixTraj   = NULL;
    _tpart       = TrkParticle::e_minus;        // 11
    _fdir        = TrkFitDirection::downstream; // = 0
    _shcol       = 0;
    _shpos       = 0;
    _shfcol      = 0;
    _krep        = 0;
    _fit         = TrkErrCode::fail;
    _nt0iter     = 0;
    _nweediter   = 0;
    _nunweediter = 0;

    _hitIndices  = new vector<StrawHitIndex>;

    static unsigned icount(0);

    while(++icount<10){
      std::cout << "KalFitResultNew is DEPRECATED and this constructor should never be called.  " 
		<< "This message will repeat " << 10-icount << " more times." << std::endl;
    }
  }


//-----------------------------------------------------------------------------
  KalFitResultNew::~KalFitResultNew() {
    if (_helixTraj) delete _helixTraj;
    delete _hitIndices;
  }

//-----------------------------------------------------------------------------
  void KalFitResultNew::deleteTrack() {
    if(_krep != 0){
      delete _krep;
      _krep = NULL; 
    }
  } 

//-----------------------------------------------------------------------------
  void KalFitResultNew::init() {
    deleteTrack();
    _hitIndices->clear();
    _missingHits.clear();
    _doca.clear();
  } 

//-----------------------------------------------------------------------------
  KalRep*  KalFitResultNew::stealTrack() { 
    KalRep* retval = _krep; 
    _krep          = 0; 
    return retval; 
  }

}

