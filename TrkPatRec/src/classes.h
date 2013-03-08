//
// Build a dictionary.
//
// $Id: classes.h,v 1.4 2013/03/08 04:33:26 brownd Exp $
// $Author: brownd $
// $Date: 2013/03/08 04:33:26 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "TrkPatRec/inc/TrkHitFilter.hh"
#include "TrkPatRec/inc/StrawHitInfo.hh"

template class std::vector<mu2e::TrkHitFilter>;
template class std::vector<mu2e::StrawHitInfo>;

