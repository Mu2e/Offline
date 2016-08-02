//
// Original author B.Echenard
//

#ifndef CaloCluster_CrystalContentMC_HH_
#define CaloCluster_CrystalContentMC_HH_


// Mu2e includes
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "MCDataProducts/inc/CaloShower.hh"
#include "MCDataProducts/inc/CaloHitMCTruthAssn.hh"
#include "CaloMC/inc/CaloContentSim.hh"



// C++ includes
#include <map>
#include <iostream>


namespace mu2e {


    class CrystalContentMC {


	 public:

             typedef std::map<art::Ptr<SimParticle>, CaloContentSim> SimContentMap;

             CrystalContentMC(const Calorimeter& cal, const CaloHitMCTruthAssns& caloHitTruth, const CaloCrystalHit& caloCrystalHit);
             ~CrystalContentMC(){};
	     	     

       	     const SimContentMap& simContentMap() const {return _simContentMap;}
	     const bool           hasConversion() const {return _hasConversion;} 
	     const double         eDepTot()       const {return _eDepTot;}
	     const double         time()          const {return _time;}
	     


	 private:

             void fillCrystal(const Calorimeter& cal, const CaloHitMCTruthAssns& caloHitTruth, const CaloCrystalHit& caloCrystalHit);
	     
	     SimContentMap _simContentMap;
	     bool          _hasConversion;
	     double        _eDepTot;
	     double        _time;

    };


} 

#endif
