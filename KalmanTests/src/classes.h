//
// Build a dictionary.
//
// $Id: classes.h,v 1.6 2014/04/18 16:40:18 kutschke Exp $
// $Author: kutschke $
// $Date: 2014/04/18 16:40:18 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "art/Persistency/Common/Wrapper.h"
#include "CLHEP/Vector/ThreeVector.h"

using namespace CLHEP;
#include "KalmanTrack/KalRep.hh"
#include "KalmanTests/inc/KalRepCollection.hh"
#include "KalmanTests/inc/KalRepPtrCollection.hh"
#include "KalmanTests/inc/KalFitMC.hh"

template class art::Wrapper<mu2e::KalRepCollection>;
template class art::Wrapper<mu2e::KalRepPtrCollection>;
template class std::vector<mu2e::TrkStrawHitInfo>;
