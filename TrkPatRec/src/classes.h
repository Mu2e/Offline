//
// Build a dictionary.
//
// $Id: classes.h,v 1.5 2014/05/31 14:28:10 brownd Exp $
// $Author: brownd $
// $Date: 2014/05/31 14:28:10 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "TrkPatRec/inc/TrkPatRec.hh"
#include "TrkPatRec/inc/TrkHitFilter.hh"
#include "TrkPatRec/inc/StrawHitInfo.hh"

template class std::vector<mu2e::TrkHitFilter>;
template class std::vector<mu2e::StrawHitInfo>;
template class std::vector<mu2e::TimePeakHitInfo>;

