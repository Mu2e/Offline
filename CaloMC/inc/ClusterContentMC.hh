//
// Original author B.Echenard
//

#ifndef CaloCluster_ClusterContentMC_HH_
#define CaloCluster_ClusterContentMC_HH_


// Mu2e includes
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "MCDataProducts/inc/CaloClusterMCTruthAssn.hh"
#include "MCDataProducts/inc/CaloShower.hh"
#include "CaloMC/inc/CaloContentSim.hh"


#include "CLHEP/Vector/ThreeVector.h"


// C++ includes
#include <map>
#include <iostream>


namespace mu2e {


    class ClusterContentMC {


	 public:

             typedef std::map<art::Ptr<SimParticle>, CaloContentSim> SimContentMap;

             ClusterContentMC(const Calorimeter& cal, const CaloClusterMCTruthAssns& caloClusterTruth, const CaloCluster& cluster);
             ~ClusterContentMC(){};
	     	     

       	     const SimContentMap& simContentMap() const {return _simContentMap;}
	     const bool           hasConversion() const {return _hasConversion;} 
	     const double         eDepTot()       const {return _eDepTot;}
	     const double         time()          const {return _time;}
	     


	 private:

             void fillCluster(const Calorimeter& cal, const CaloClusterMCTruthAssns& caloClusterTruth, const CaloCluster& cluster);

	     SimContentMap _simContentMap;
	     bool          _hasConversion;
	     double        _eDepTot;
	     double        _time;

    };


} 

#endif
