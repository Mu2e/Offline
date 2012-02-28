//
//
//
// $Id: CaloClusterFinder.hh,v 1.1 2012/02/28 22:24:48 gianipez Exp $
// $Author: gianipez $
// $Date: 2012/02/28 22:24:48 $
//
// Original author G. Pezzullo & G. Tassielli
//


#ifndef CALOCLUSTERFINDER_HH_
#define CALOCLUSTERFINDER_HH_

#include "CaloCluster/inc/CaloClusterUtilities.hh"
#include "RecoDataProducts/inc/CaloCrystalHit.hh"
//#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"

//#include "RecoDataProducts/inc/CaloCluster.hh"
//#include "RecoDataProducts/inc/CaloClusterCollection.hh"

namespace mu2e {

//define a map which key is the index of the row "R" of the vane, and the contained object is an other map which key is the column index "Z" and also contain a vector.
//The vector contains pairs of "CaloCrystalHit" and a index which is the position of the "CaloCrystalHit" in the vector "CaloCrystalHitCollection" generated in the event
typedef std::map<unsigned int, std::map<unsigned int, std::vector<std::pair<CaloCrystalHit *, size_t > > >  > MatrixCaloHit;

//define a map which key is the vane's index and contains object of type "CaloCrystalHit". In that way we have a complete topology description of the calorimeter
typedef std::map<unsigned int, MatrixCaloHit> VanesMap;


//define the object in which we store a single cluster. the key is the vane's index and the pair contains (row, column)
typedef std::multimap<unsigned int, std::pair<unsigned int, unsigned int> > ClusterData;//row, cloumn, hitId


class CaloClusterFinder {
public :
        CaloClusterFinder(CaloClusterParameters *ccl):
        HitTimeMin(0.0)
        {
                _ccl = ccl;
                looked = new bool *[(unsigned int) _ccl->_nRow];
                for(int i=0; i<  _ccl->_nRow; ++i){
                        looked[i]= new bool [(unsigned int) _ccl->_nColum];
                }

        }
        ~CaloClusterFinder (){
                for(int i=0; i< _ccl->_nRow; ++i){
                        delete [] looked[i];
                }
                delete [] looked;
        }

        virtual bool find(ClusterData &cluster, MatrixCaloHit &data , unsigned int row, unsigned int colum ){
                return false;
        }

        void setFirstHitTime(float &time){
                HitTimeMin = time;
        }

        virtual void initializeNewSearch(){
                for(int y=0; y<_ccl->_nRow; ++y){
                        for(int u=0; u<_ccl->_nColum; ++u){
                                looked[y][u]=false;
                        }
                }
        }

protected:
        CaloClusterParameters *_ccl;
        float HitTimeMin;
        bool **looked;
};

}



#endif /* CALOCLUSTERFINDER_HH_ */
