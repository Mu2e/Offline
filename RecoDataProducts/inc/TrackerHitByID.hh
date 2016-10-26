#ifndef RecoDataProducts_TrackerHitByID_hh
#define RecoDataProducts_TrackerHitByID_hh
//
// map of the Tracker hit to be accessed by ID
//
// $Id: TrackerHitByID.hh,v 1.1 2012/04/20 07:18:34 tassiell Exp $
// $Author: tassiell $
// $Date: 2012/04/20 07:18:34 $
//
// Original author G. Tassielli
//

// C++ includes
#include <map>
//#include <algorithm>

// Mu2e includes
#include "RecoDataProducts/inc/StrawHit.hh"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

  typedef std::multimap<unsigned long int, art::Ptr<mu2e::StrawHit>, std::less<unsigned long int> > TrackerHitByID;

} // namespace mu2e

#endif /* RecoDataProducts_TrackerHitByID_hh */
