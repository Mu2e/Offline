#ifndef RecoDataProducts_StrawHitFlagCollection_hh
#define RecoDataProducts_StrawHitFlagCollection_hh

//
// Define a type for a collection of StrawHitFlag objects.
//
// $Id: StrawHitFlagCollection.hh,v 1.1 2013/03/08 04:29:49 brownd Exp $
// $Author: brownd $
// $Date: 2013/03/08 04:29:49 $
//
// Original author David Brown
//

#include <vector>

#include "RecoDataProducts/inc/StrawHitFlag.hh"

namespace mu2e {
  typedef std::vector<mu2e::StrawHitFlag> StrawHitFlagCollection;
}

#endif /* RecoDataProducts_StrawHitFlagCollection_hh */
