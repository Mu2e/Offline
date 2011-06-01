//
// StrawCluster to be created based on CaloHit's
//
// $Id: StrawCluster.cc,v 1.1 2011/06/01 21:13:31 wenzel Exp $
// $Author: wenzel $
// $Date: 2011/06/01 21:13:31 $
//

// C++ includes
#include <ostream>

// Framework includes.
//#include "art/Utilities/Exception.h"
#include "art/Persistency/Provenance/ProductID.h"

// Mu2e includes
#include "RecoDataProducts/inc/StrawCluster.hh"

using namespace std;

namespace mu2e {

  StrawCluster::StrawCluster(std::vector<DPIndex>&  hitIndices)
  {
     _StrawHitIndices=hitIndices;
  }
} // namespace mu2e
