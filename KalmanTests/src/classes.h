//
// Build a dictionary.
//
// $Id: classes.h,v 1.7 2014/04/26 18:13:38 murat Exp $
// $Author: murat $
// $Date: 2014/04/26 18:13:38 $
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

template class art::Ptr<KalRep>;
template class art::Wrapper<mu2e::KalRepCollection>;
template class art::Wrapper<mu2e::KalRepPtrCollection>;
template class std::vector<mu2e::TrkStrawHitInfo>;
