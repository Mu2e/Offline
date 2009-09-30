#ifndef ToyDP_ToyHitCollection_hh
#define ToyDP_ToyHitCollection_hh

//
// Define a type for a collection of ToyHits.
//
// $Id: ToyHitCollection.hh,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "ToyDP/inc/ToyHit.hh"

namespace mu2e {
   typedef std::vector<mu2e::ToyHit> ToyHitCollection;
}

#endif
