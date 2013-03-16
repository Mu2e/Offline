//
// Build a dictionary.
//
// $Id: classes.h,v 1.5 2013/03/16 04:27:45 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/16 04:27:45 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "art/Persistency/Common/Wrapper.h"
#include "CLHEP/Vector/ThreeVector.h"

using namespace CLHEP;
#include "KalmanTrack/KalRep.hh"
#include "KalmanTests/inc/KalRepCollection.hh"
#include "KalmanTests/inc/KalFitMC.hh"

template class art::Wrapper<mu2e::KalRepCollection>;
template class std::vector<mu2e::TrkStrawHitInfo>;
