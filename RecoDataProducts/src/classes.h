//
// Build a dictionary.
//
// $Id: classes.h,v 1.3 2011/06/07 21:32:22 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/06/07 21:32:22 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "art/Persistency/Common/Wrapper.h"

#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawClusterCollection.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/HoughCircleCollection.hh"
#include "RecoDataProducts/inc/SubEventCollection.hh"

// Cannot use the typedefs in here - not sure why.
template class art::Ptr<mu2e::CaloHit>;
template class std::vector<art::Ptr<mu2e::CaloHit> >;
template class art::Ptr<mu2e::StrawHit>;
template class std::vector<art::Ptr<mu2e::StrawHit> >;

template class art::Wrapper<mu2e::StrawHitCollection>;
template class art::Wrapper<mu2e::StrawClusterCollection>;
template class art::Wrapper<mu2e::CaloHitCollection>;
template class art::Wrapper<mu2e::CaloCrystalHitCollection>;
template class art::Wrapper<mu2e::HoughCircleCollection>;
template class art::Wrapper<mu2e::SubEventCollection>;
