#ifndef CaloCluster_ClusterContentMC_HH_
#define CaloCluster_ClusterContentMC_HH_


#include "CaloMC/inc/CaloContentSim.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/CaloClusterMCTruthAssn.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"

#include <map>
#include <iostream>


namespace mu2e {


    class ClusterContentMC {


	 public:

             typedef std::map<art::Ptr<SimParticle>, CaloContentSim> SimContentMap;

             ClusterContentMC(const Calorimeter& cal, const CaloClusterMCTruthAssns& caloClusterTruth, const CaloCluster& cluster);
             ~ClusterContentMC(){};
	     	     

       	     const SimContentMap& simContentMap() const {return simContentMap_;}
	     const bool           hasConversion() const {return hasConversion_;} 
	     const double         eDepTot()       const {return eDepTot_;}
	     const double         time()          const {return time_;}
	     


	 private:

             void fillCluster(const Calorimeter& cal, const CaloClusterMCTruthAssns& caloClusterTruth, const CaloCluster& cluster);

	     SimContentMap simContentMap_;
	     bool          hasConversion_;
	     double        eDepTot_;
	     double        time_;

    };


} 

#endif
