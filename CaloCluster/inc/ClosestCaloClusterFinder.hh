//
//
//
// $Id: ClosestCaloClusterFinder.hh,v 1.1 2012/02/28 22:24:48 gianipez Exp $
// $Author: gianipez $
// $Date: 2012/02/28 22:24:48 $
//
// Original author G. Pezzullo & G. Tassielli
//


#ifndef CLOSESTCALOCLUSTERFINDER_HH_
#define CLOSESTCALOCLUSTERFINDER_HH_

#include "CaloCluster/inc/CaloClusterFinder.hh"

namespace mu2e {

class ClosestCaloClusterFinder : public CaloClusterFinder {
public :
        ClosestCaloClusterFinder(CaloClusterParameters *ccl):
        CaloClusterFinder::CaloClusterFinder(ccl){
        }
        ~ClosestCaloClusterFinder (){}

        bool find(ClusterData &cluster, MatrixCaloHit &data , unsigned int row, unsigned int colum );
};
}


#endif /* CLOSESTCALOCLUSTERFINDER_HH_ */
