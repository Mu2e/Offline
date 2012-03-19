//
// General utilities for the calorimeter's studies
//
// $Id: CaloClusterUtilities.hh,v 1.3 2012/03/19 19:35:42 gianipez Exp $
// $Author: gianipez $
// $Date: 2012/03/19 19:35:42 $
//
// Original author G. Pezzullo & G. Tassielli & G. Onorato
//

//art includes
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Principal/Handle.h"

//Mu2e includes
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
//#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHit.hh"
#include "CLHEP/Vector/ThreeVector.h"

//---------------
// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Utilities/Exception.h"
#include "art/Framework/Principal/Event.h"

// Mu2e includes

#include "MCDataProducts/inc/SimParticleCollection.hh"

//------

//c++ includes
#include <string>

//Root includes
#include "TMath.h"


using namespace std;

namespace mu2e {


class MCCaloUtilities {

public:

  MCCaloUtilities();

  ~MCCaloUtilities();

  void setTrackAndRO(const art::Event & event,
                     std::string const &_g4ModuleLabel,
                     SimParticleCollection::key_type track,
                     unsigned RO);

  void printOutCaloInfo();

  bool fromOutside();

  bool primary();

  bool generated();

  int startingVane();

  int getStartingVane(CLHEP::Hep3Vector origin);

  int localVane();

private:

  unsigned _localRO;
  unsigned _localCrystal;
  unsigned _localVane;
  int _startingVane;
  bool _fromOutside, _primary, _generated;

};

//--------------------------------------------------------------------------------------------------------------

std::string & TOUpper(std::string &in);

double cry(double val);

//define a map which key is the index of the row "R" of the vane, and the contained object is an other map which key is the column index "Z" and also contain a vector.
//The vector contains pairs of "CaloCrystalHit" and a index which is the position of the "CaloCrystalHit" in the vector "CaloCrystalHitCollection" generated in the event
typedef std::map<unsigned int, std::map<unsigned int, std::vector<std::pair<CaloCrystalHit *, size_t > > >  > MatrixCaloHit;

//define a map which key is the vane's index and contains object of type "CaloCrystalHit". In that way we have a complete topology description of the calorimeter
typedef std::map<unsigned int, MatrixCaloHit> VanesMap;


//define the object in which we store a single cluster. the key is the row index, and the pair contains the column index and the position of the CaloCrystalHit in the vector
//stored in the container "vanesMap"
typedef std::multimap<unsigned int, std::pair<unsigned int, unsigned int> > ClusterData;//row, cloumn, hitId


struct ClusterMap{
        int _cluSize;
        int _COGcrySize;
        std::vector<double> _rowVec;
        std::vector<double> _columnVec;
        std::vector<double> _cryEdepVec;
        std::vector<double> _COGrowVec;
        std::vector<double> _COGcolumnVec;
        int _vane;
        int _cluCogRow;
        int _cluCogColumn;
        CLHEP::Hep3Vector _cluCOG;
        float   _showerDir;
        float   _errShowerDir;


//        bool operator<( const elecData other) const{
//                return ( _impTime< other._impTime);
//        }
        ClusterMap & operator=(const ClusterMap& other) {
                _cluSize        = other._cluSize;
                _COGcrySize     = other._COGcrySize;
                _cluSize        = other._cluSize;
                _vane           = other._vane;
                _cluCogRow      = other._cluCogRow;
                _cluCogColumn   = other._cluCogColumn;
                _cluCOG         = other._cluCOG;
                _rowVec         = other._rowVec;
                _columnVec      = other._columnVec;
                _cryEdepVec     = other._cryEdepVec;
                _COGrowVec      = other._COGcolumnVec;
                _COGcolumnVec   = other._COGcolumnVec;
                _showerDir      = other._showerDir;
//                _cluCOG.x() = other._cluCOG.x();
//                _cluCOG.y()    = other._cluCOG.y();
//                _cluCOG.z() = other._cluCOG.z();
//                _rowVec.size() = other._rowVec.size();
//                _columnVec.size() = other._columnVec.size();
//                _COGrowVec.size() = other._COGcolumnVec.size();
//                _COGcolumnVec.size()    = other._COGcolumnVec.size();



                return *this;
        }
        ClusterMap():
                _cluSize(0),
                _cluCogRow(0.0),
                _cluCogColumn(0.0){
        }
};


//parameters used to build the cluster
class CaloClusterParameters{
public:
        CaloClusterParameters(){}

        double                                _deltaTime;//[ns] time window requested to the crystals of each cluster
        int                             _nCryPerCluster;// minimum number of crystals for defining a cluster
        double                                _EnoiseCut;//[MeV] equal to 3 sigma noise
        double                              _EclusterCut;
        int                                  _nRow;
        int                                _nColum ;

        ~CaloClusterParameters(){}

};

//the following procedure returns the point of impact (3Vector) calculated using the crystal's fraction energy as the weight of the mean
void cog(CaloCluster &cluster);

void cog_depth(CaloCluster &cluster, double depth, ClusterMap &clusterMap);

//on the following we implement an algorithm for the cog which uses the logarithm of the energy as weight (w_{i}). This is not correct, as references show, we need to calculate from simulation an offset to add at each w_{i}
//void LOGcog(CaloCluster &cluster);
CLHEP::Hep3Vector LOGcog(CaloCluster &cluster, double w, double depth);

void LOGcogMap(CaloCluster &cluster, double w, double depth, ClusterMap &clusterMap );

}

//#endif /* CaloClusterUtilities_hh */
