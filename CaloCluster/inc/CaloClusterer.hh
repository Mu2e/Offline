//
// Functions used to build the clusters
//
// $Id: CaloClusterer.hh,v 1.3 2013/03/15 15:52:03 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/15 15:52:03 $
//
// Original author G. Pezzullo & G. Tassielli
//



#ifndef CALOCLUSTERER_HH_
#define CALOCLUSTERER_HH_

// C++ includes
#include <string>
#include "CaloCluster/inc/ClosestClusterFinder.hh"
#include "RecoDataProducts/inc/CaloCrystalHit.hh"


namespace mu2e {

//define a map which key is the index of the row "R" of the vane, and the contained object is an other map which key is the column index "Z" and also contain a vector.
//The vector contains pairs of "CaloCrystalHit" and a index which is the position of the "CaloCrystalHit" in the vector "CaloCrystalHitCollection" generated in the event
typedef std::map<unsigned int, std::map<unsigned int, std::vector<std::pair<CaloCrystalHit *, size_t > > >  > MatrixCaloHit;

//define a map which key is the vane's index and contains object of type "CaloCrystalHit". In that way we have a complete topology description of the calorimeter
typedef std::map<unsigned int, MatrixCaloHit> VanesMap;


//define the object in which we store a single cluster. the key is the vane's index and the pair contains (row, column)
typedef std::multimap<unsigned int, std::pair<unsigned int, unsigned int> > ClusterData;//row, cloumn, hitId


class CaloClusterer{

private:

        std::unique_ptr<ClusterFinder> _ccf;
        CaloClusterParameters *_param;

public:

        CaloClusterer(int                                  nRow = 0,
                        int                               nColum = 0,
                        std::string flagAlgo = "") {
                _param = new CaloClusterParameters();
                _param->_nRow = nRow;
                _param->_nColum = nColum;
                _param->_deltaTime = 0.0;
                _param->_nCryPerCluster = 0;
                _param->_EnoiseCut = 0.0;
                _param->_EclusterCut = 0.0;

                if(flagAlgo.compare("CLOSEST") == 0){
                        _ccf.reset( new ClosestClusterFinder(_param) );
                }else if(flagAlgo.compare("OTHERalgorithm")){
                        //_ccf = new OtherAlgorithmClusterFinder(_param);
                }else {
                        _ccf.reset( new ClusterFinder(_param) );
                }

        }

        CaloClusterer(int                                  nRow,
                        int                               nColum,
                        double                               deltaTime,
                        int                           nCryPerCluster,
                        double                             EnoiseCut,
                        double                           EclusterCut,
                        std::string                      flagAlgo=""){
                _param = new CaloClusterParameters();
                _param->_nRow = nRow;
                _param->_nColum = nColum;
                _param->_deltaTime = deltaTime;
                _param->_nCryPerCluster = nCryPerCluster;
                _param->_EnoiseCut = EnoiseCut;
                _param->_EclusterCut = EclusterCut;

                if(flagAlgo.compare("CLOSEST") == 0){
                        _ccf.reset( new ClosestClusterFinder(_param) );
                }else if(flagAlgo.compare("OTHERalgorithm")){
                        //_ccf = new OtherAlgorithmClusterFinder(_param);
                }else {
                        _ccf.reset( new ClusterFinder(_param) );
                }
        }

        bool find(ClusterData &cluster, MatrixCaloHit &data , unsigned int row, unsigned int colum,int hitid=-1 ){
                return _ccf->find(cluster, data , row, colum, hitid);
        }

        void setFirstHitTime(float &time){
                _ccf->setFirstHitTime(time);
        }

        void setFirstHitTime(float time){
                _ccf->setFirstHitTime(time);
        }

        void initializeNewSearch(){
                _ccf->initializeNewSearch();
        }

};

} // end namespace mu2e

#endif
