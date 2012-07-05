//
// Build a dictionary.
//
// $Id: classes.h,v 1.2 2012/07/05 21:38:53 brownd Exp $
// $Author: brownd $
// $Date: 2012/07/05 21:38:53 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "TrkPatRec/inc/TrkHitFilter.hh"
#include "TrkPatRec/inc/StrawHitInfo.hh"

template class std::vector<mu2e::TrkHitFilter>;
template class std::vector<mu2e::StrawHitInfo>;

