//
// Definition of one algorithm used to build the clusters. It's based on the proximity principal
//
// $Id: ClosestCaloClusterFinder.cc,v 1.2 2012/03/07 18:00:38 gianipez Exp $
// $Author: gianipez $
// $Date: 2012/03/07 18:00:38 $
//
// Original author G. Pezzullo & G. Tassielli
//

#include "CaloCluster/inc/ClosestCaloClusterFinder.hh"

namespace mu2e{

bool ClosestCaloClusterFinder::find( ClusterData &cluster, MatrixCaloHit &data , unsigned int row, unsigned int colum, int hitid ){


        //physical parameters used for the clustering operation
        double deltaTime = _ccl->_deltaTime;//[ns] time window requested to the crystals of each cluster
        int nCryPerCluster = _ccl->_nCryPerCluster;// minimum number of crystals for defining a cluster
        double noiseCut = _ccl->_EnoiseCut;//[MeV] equal to 3 sigma noise


        unsigned int tmp_row, tmp_colum;
        tmp_row=row;
        tmp_colum=colum;

        if (hitid!=-1) {
                cluster.insert(ClusterData::value_type( tmp_row, std::pair<unsigned int, unsigned int>( tmp_colum, hitid ) ));
                looked[tmp_row][tmp_colum] = true;
        }

        for (int irw=-1; irw<2; irw++) {
                tmp_row=row+irw;
                MatrixCaloHit::iterator ida = data.find(tmp_row);
                if ( ida != data.end()  ) {
                        for (int icl=-1; icl<2; icl++) {
                                tmp_colum=colum+icl;
                                std::map<unsigned int, std::vector<std::pair<CaloCrystalHit *, size_t > > >::iterator ida2 = ida->second.find(tmp_colum);
                                if ( (ida2 != ida->second.end()) && !looked[tmp_row][tmp_colum] ) {

                                        looked[tmp_row][tmp_colum]=true;

                                        for(std::vector<std::pair<CaloCrystalHit *, size_t > >::iterator jHit = ida2->second.begin(); jHit != ida2->second.end(); ++jHit ){
                                                if(fabs(HitTimeMin - jHit->first->time()) <= deltaTime  && fabs(HitTimeMax - jHit->first->time()) <= deltaTime  && jHit->first->energyDep() >= noiseCut ){
                                                        if(HitTimeMin > jHit->first->time()) HitTimeMin = jHit->first->time();
                                                        if(HitTimeMax < jHit->first->time()) HitTimeMax = jHit->first->time();
                                                        cluster.insert( ClusterData::value_type( tmp_row, std::pair<unsigned int, unsigned int>( tmp_colum, jHit-ida2->second.begin()/*jHit->second*/) ));
                                                        find(cluster, data, tmp_row, tmp_colum);
                                                        break;
                                                }
                                        }
                                }
                        }//end for(icl)
                }
        }//end for(irw)

        return (cluster.size() >= (size_t) nCryPerCluster);
}

}
