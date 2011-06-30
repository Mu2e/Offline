//
// Build a dictionary.
//
// $Id: classes.h,v 1.3 2011/06/30 20:46:49 mu2ecvs Exp $
// $Author: mu2ecvs $
// $Date: 2011/06/30 20:46:49 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "art/Persistency/Common/Wrapper.h"

using namespace CLHEP;
#include "TrkBase/TrkRecoTrk.hh"
#include "KalmanTests/inc/TrkRecoTrkCollection.hh" 
#include "KalmanTests/inc/KalFitMC.hh"

template class art::Wrapper<mu2e::TrkRecoTrkCollection>;
template class std::vector<mu2e::TrkStrawHitInfo>;
