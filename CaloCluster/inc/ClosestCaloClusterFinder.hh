//
//
//
// $Id: ClosestCaloClusterFinder.hh,v 1.2 2012/03/07 18:00:38 gianipez Exp $
// $Author: gianipez $
// $Date: 2012/03/07 18:00:38 $
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

        bool find(ClusterData &cluster, MatrixCaloHit &data , unsigned int row, unsigned int colum, int hitid=-1 );
};
}


#endif /* CLOSESTCALOCLUSTERFINDER_HH_ */
