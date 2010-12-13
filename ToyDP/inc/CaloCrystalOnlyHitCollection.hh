#ifndef ToyDP_CaloCrystalOnlyHitCollection_hh
#define ToyDP_CaloCrystalOnlyHitCollection_hh

//
// Define a type for a collection of CaloCrystalOnlyHit objects.
//
// $Id: CaloCrystalOnlyHitCollection.hh,v 1.1 2010/12/13 06:10:33 logash Exp $
// $Author: logash $
// $Date: 2010/12/13 06:10:33 $
//
// Original author Ivan Logashenko
//

#include <vector>

#include "ToyDP/inc/CaloCrystalOnlyHit.hh"

namespace mu2e {
   typedef std::vector<mu2e::CaloCrystalOnlyHit> CaloCrystalOnlyHitCollection;
}

#endif
