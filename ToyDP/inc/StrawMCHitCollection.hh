#ifndef ToyDP_StrawMCHitCollection_hh
#define ToyDP_StrawMCHitCollection_hh

//
// Define a type for a collection of StrawMCHits.
//
// $Id: StrawMCHitCollection.hh,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "ToyDP/inc/StrawMCHit.hh"

namespace mu2e {
   typedef std::vector<mu2e::StrawMCHit> StrawMCHitCollection;
}

#endif
