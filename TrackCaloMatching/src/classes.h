//
// Build a dictionary.
//
// $Id: classes.h,v 1.1 2012/07/10 00:02:19 gianipez Exp $
// $Author: gianipez $
// $Date: 2012/07/10 00:02:19 $
//
// Original author G. Pezzullo
//

#include <vector>

#include "art/Persistency/Common/Wrapper.h"

#include "TrackCaloMatching/inc/TrkToCaloExtrapol.hh"
#include "TrkBase/TrkRecoTrk.hh"
#include "TrackCaloMatching/inc/TrkToCaloExtrapolCollection.hh"
//#include "KalmanTests/inc/KalFitMC.hh"

template class art::Ptr<const TrkRecoTrk * const >;

template class art::Ptr<mu2e::TrkToCaloExtrapol>;
template class std::vector<art::Ptr<mu2e::TrkToCaloExtrapol> >;

template class art::Wrapper<mu2e::TrkToCaloExtrapolCollection>;
