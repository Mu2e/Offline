#ifndef ToyDP_CrudeStrawHitPData_hh
#define ToyDP_CrudeStrawHitPData_hh

//
// Define the persistent data for a collection of 
// CrudeStrawHit objects.
//
// $Id: CrudeStrawHitPData.hh,v 1.1 2009/10/22 15:51:23 kutschke Exp $
// $Author: kutschke $
// $Date: 2009/10/22 15:51:23 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "ToyDP/inc/CrudeStrawHit.hh"

namespace mu2e {
  typedef std::vector<mu2e::CrudeStrawHit> CrudeStrawHitPData;
}

#endif
