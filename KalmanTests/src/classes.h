//
// Build a dictionary.
//
// $Id: classes.h,v 1.10 2014/08/22 20:51:05 brownd Exp $
// $Author: brownd $
// $Date: 2014/08/22 20:51:05 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "art/Persistency/Common/Wrapper.h"
#include "CLHEP/Vector/ThreeVector.h"

using namespace CLHEP;
#include "BTrk/KalmanTrack/KalRep.hh"
#include "KalmanTests/inc/KalRepCollection.hh"
#include "KalmanTests/inc/TrkStrawHitInfo.hh"
#include "KalmanTests/inc/KalFitResult.hh"
#include "KalmanTests/inc/PanelAmbigStructs.hh"

template class art::Wrapper<mu2e::KalRepCollection>;
template class std::vector<mu2e::TrkStrawHitInfo>;
template class std::vector<mu2e::TrkStrawHitInfoMC>;
template class std::vector<mu2e::TrkStrawHitInfo_old>;
template class std::vector<mu2e::PanelAmbig::PanelResult>;
template class std::vector<mu2e::PanelAmbig::TSHUInfo>;
template class std::vector<mu2e::PanelAmbig::TSHMCUInfo>;

template class art::Wrapper<mu2e::KalFitResultCollection>;
