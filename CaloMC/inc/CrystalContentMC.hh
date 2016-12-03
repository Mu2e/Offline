#ifndef CaloCluster_CrystalContentMC_HH_
#define CaloCluster_CrystalContentMC_HH_


#include "CaloMC/inc/CaloContentSim.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/CaloHitMCTruthAssn.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"

#include <map>
#include <iostream>


namespace mu2e {


    class CrystalContentMC {


	 public:

             typedef std::map<art::Ptr<SimParticle>, CaloContentSim> SimContentMap;

             CrystalContentMC(const Calorimeter& cal, const CaloHitMCTruthAssns& caloHitTruth, const CaloCrystalHit& caloCrystalHit);
             ~CrystalContentMC(){};
	     	     

       	     const SimContentMap& simContentMap() const {return simContentMap_;}
	     const bool           hasConversion() const {return hasConversion_;} 
	     const double         eDepTot()       const {return eDepTot_;}
	     const double         time()          const {return time_;}
	     


	 private:

             void fillCrystal(const Calorimeter& cal, const CaloHitMCTruthAssns& caloHitTruth, const CaloCrystalHit& caloCrystalHit);
	     
	     SimContentMap simContentMap_;
	     bool          hasConversion_;
	     double        eDepTot_;
	     double        time_;

    };


} 

#endif
