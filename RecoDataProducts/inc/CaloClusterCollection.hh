//
// Calorimeter's data clusters container
//
// $Id: CaloClusterCollection.hh,v 1.1 2012/02/28 22:26:01 gianipez Exp $
// $Author: gianipez $
// $Date: 2012/02/28 22:26:01 $
//
// Original author G. Pezzullo & G. Tassielli
//


#ifndef RecoDataProducts_CaloClusterCollection_hh
#define RecoDataProducts_CaloClusterCollection_hh



#include <vector>

#include "RecoDataProducts/inc/CaloCluster.hh"

namespace mu2e {
  typedef std::vector<mu2e::CaloCluster> CaloClusterCollection;
}

#endif /* RecoDataProducts_CaloClusterCollection_hh */
