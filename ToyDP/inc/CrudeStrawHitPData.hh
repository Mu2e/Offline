#ifndef ToyDP_CrudeStrawHitPData_hh
#define ToyDP_CrudeStrawHitPData_hh

//
// Define the persistent data for a collection of 
// CrudeStrawHit objects.
//
// $Id: CrudeStrawHitPData.hh,v 1.2 2011/05/17 15:41:36 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:36 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "ToyDP/inc/CrudeStrawHit.hh"

namespace mu2e {
  typedef std::vector<mu2e::CrudeStrawHit> CrudeStrawHitPData;
}

#endif /* ToyDP_CrudeStrawHitPData_hh */
