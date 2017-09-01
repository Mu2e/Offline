//
// Struct to hold BaBar Kalman fit
//
// $Id: KalFitResult.cc,v 1.1 2012/08/31 23:21:02 brownd Exp $
// $Author: brownd $ 
// $Date: 2012/08/31 23:21:02 $
//

// the following has to come before other BaBar includes
#include "BTrk/BaBar/BbrCollectionUtils.hh"
#include "CalPatRec/inc/KalFitResult.hh"
namespace mu2e {

  KalFitResult::KalFitResult() {
    _tdef        = 0;
    _shcol       = 0;
    _krep        = 0;
    _fit         = TrkErrCode::fail;
    _nt0iter     = 0;
    _nweediter   = 0;
    _nunweediter = 0;
    _ninter      = 0;
    static unsigned icount(0);
    while(++icount<10){
      std::cout << "KalFitResult is DEPRECATED and this constructor should never be called.  " 
		<< "This message will repeat " << 10-icount << " more times." << std::endl;
    }
  }

  KalFitResult::KalFitResult(const TrkDefHack* hdef) {
    _tdef        = hdef;
    _shcol       = hdef->shcol();
    _krep        = 0;
    _fit         = TrkErrCode::fail;
    _nt0iter     = 0;
    _nweediter   = 0;
    _nunweediter = 0;
    _ninter      = 0;
      
    static unsigned icount(0);
    while(++icount<10){
      std::cout << "KalFitResult is DEPRECATED: the code calling this function needs to be refactored.  " 
		<< "This message will repeat " << 10-icount << " more times." << std::endl;
    }
  }
    
  KalFitResult::~KalFitResult() {
  }

//-----------------------------------------------------------------------------
// 
  void KalFitResult::deleteTrack() {
    if(_krep != 0){
      _hits.clear();
      delete _krep;
      _krep = NULL; 
    }
    else {
      std::for_each(_hits.begin(),_hits.end(),babar::Collection::DeleteObject());
    }
  } 
}

