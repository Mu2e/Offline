//
// Build a dictionary.
//
// $Id: classes.h,v 1.4 2012/07/23 17:52:27 brownd Exp $
// $Author: brownd $
// $Date: 2012/07/23 17:52:27 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "art/Persistency/Common/Wrapper.h"

using namespace CLHEP;
#include "KalmanTrack/KalRep.hh"
#include "KalmanTests/inc/KalRepCollection.hh" 
#include "KalmanTests/inc/KalFitMC.hh"

template class art::Wrapper<mu2e::KalRepCollection>;
template class std::vector<mu2e::TrkStrawHitInfo>;
