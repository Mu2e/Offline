#ifndef ToyDP_CrudeStrawHitPData_hh
#define ToyDP_CrudeStrawHitPData_hh

//
// Define the persistent data for a collection of
// CrudeStrawHit objects.
//
// $Id: CrudeStrawHitPData.hh,v 1.3 2011/05/18 02:27:19 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:19 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "ToyDP/inc/CrudeStrawHit.hh"

namespace mu2e {
  typedef std::vector<mu2e::CrudeStrawHit> CrudeStrawHitPData;
}

#endif /* ToyDP_CrudeStrawHitPData_hh */
