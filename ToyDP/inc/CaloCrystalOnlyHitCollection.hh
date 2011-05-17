#ifndef ToyDP_CaloCrystalOnlyHitCollection_hh
#define ToyDP_CaloCrystalOnlyHitCollection_hh

//
// Define a type for a collection of CaloCrystalOnlyHit objects.
//
// $Id: CaloCrystalOnlyHitCollection.hh,v 1.2 2011/05/17 15:41:36 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:36 $
//
// Original author Ivan Logashenko
//

#include <vector>

#include "ToyDP/inc/CaloCrystalOnlyHit.hh"

namespace mu2e {
   typedef std::vector<mu2e::CaloCrystalOnlyHit> CaloCrystalOnlyHitCollection;
}

#endif /* ToyDP_CaloCrystalOnlyHitCollection_hh */
