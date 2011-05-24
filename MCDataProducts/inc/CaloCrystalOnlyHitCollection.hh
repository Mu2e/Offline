#ifndef MCDataProducts_CaloCrystalOnlyHitCollection_hh
#define MCDataProducts_CaloCrystalOnlyHitCollection_hh

//
// Define a type for a collection of CaloCrystalOnlyHit objects.
//
// $Id: CaloCrystalOnlyHitCollection.hh,v 1.1 2011/05/24 17:16:43 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/05/24 17:16:43 $
//
// Original author Ivan Logashenko
//

#include <vector>

#include "MCDataProducts/inc/CaloCrystalOnlyHit.hh"

namespace mu2e {
   typedef std::vector<mu2e::CaloCrystalOnlyHit> CaloCrystalOnlyHitCollection;
}

#endif /* MCDataProducts_CaloCrystalOnlyHitCollection_hh */
