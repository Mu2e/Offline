//
// Struct to hold BaBar Kalman fit
//
// $Id: KalFitResult.cc,v 1.1 2012/08/31 23:21:02 brownd Exp $
// $Author: brownd $ 
// $Date: 2012/08/31 23:21:02 $
//

// the following has to come before other BaBar includes
#include "BTrk/BaBar/BaBar.hh"
#include "TrkReco/inc/KalFitResult.hh"
namespace mu2e 
{

  void KalFitResult::deleteTrack() {
    if(_krep != 0){
      _hits.clear(); delete _krep; _krep = 0; 
//      std::cout << "deleting fit with track " << std::endl;
    } else {
// if there's no fit, we need to delete the hits
//      std::cout << "deleting " << _hits.size() << " hits from fit without track " << std::endl;
      std::for_each(_hits.begin(),_hits.end(),babar::Collection::DeleteObject());
    }
  } 
}

