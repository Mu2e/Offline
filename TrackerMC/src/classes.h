//
// Build a dictionary.
//
// $Id: classes.h,v 1.10 2014/08/22 20:51:05 brownd Exp $
// $Author: brownd $
// $Date: 2014/08/22 20:51:05 $
//

#include <vector>
#include "CLHEP/Vector/ThreeVector.h"
#include "TrackerMC/inc/SHInfo.hh"
#include "TrkChargeReco/inc/PeakFitParams.hh"

template class std::vector<mu2e::SHID>;
template class std::vector<mu2e::SHMCInfo>;
template class std::vector<mu2e::TrkChargeReco::PeakFitParams>;

