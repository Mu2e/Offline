//
// Build a dictionary.
//
// $Id: classes.h,v 1.10 2014/08/22 20:51:05 brownd Exp $
// $Author: brownd $
// $Date: 2014/08/22 20:51:05 $
//
// Original author Rob Kutschke
//

#include "art/Persistency/Common/Wrapper.h"

#include "BTrk/KalmanTrack/KalRep.hh"
#include "RecoDataProducts/inc/KalRepCollection.hh"
#include "TrkReco/inc/PanelAmbigStructs.hh"

template class art::Wrapper<mu2e::KalRepCollection>;
template class std::vector<mu2e::PanelAmbig::PanelResult>;
template class std::vector<mu2e::PanelAmbig::TSHUInfo>;
template class std::vector<mu2e::PanelAmbig::TSHMCUInfo>;

