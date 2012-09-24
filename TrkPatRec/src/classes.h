//
// Build a dictionary.
//
// $Id: classes.h,v 1.3 2012/09/24 18:39:55 brownd Exp $
// $Author: brownd $
// $Date: 2012/09/24 18:39:55 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "TrkPatRec/inc/TrkHitFilter.hh"
#include "TrkPatRec/inc/StrawHitInfo.hh"
#include "TrkPatRec/inc/DeltaHitInfo.hh"

template class std::vector<mu2e::TrkHitFilter>;
template class std::vector<mu2e::StrawHitInfo>;

