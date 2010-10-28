#ifndef ToyDP_CaloCrystalHitMCTruthCollection_hh
#define ToyDP_CaloCrystalHitMCTruthCollection_hh

//
// Define a type for a collection of CaloCrystalHitMCTruth objects.
//
// $Id: CaloCrystalHitMCTruthCollection.hh,v 1.1 2010/10/28 20:43:58 genser Exp $
// $Author: genser $
// $Date: 2010/10/28 20:43:58 $
//
// Original author Ivan Logashenko
//

#include <vector>

#include "ToyDP/inc/CaloCrystalHitMCTruth.hh"

namespace mu2e {
   typedef std::vector<mu2e::CaloCrystalHitMCTruth> CaloCrystalHitMCTruthCollection;
}

#endif
