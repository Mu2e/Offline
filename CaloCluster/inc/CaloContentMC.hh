//
// Original author B.Echenard
//

#ifndef CaloCluster_CaloContentMC_HH_
#define CaloCluster_CaloContentMC_HH_


// Mu2e includes
#include "CaloCluster/inc/CaloContentMC.hh"
#include "HitMakers/inc/CaloCrystalMCUtil.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/CaloHitSimPartMCCollection.hh"
#include "Mu2eUtilities/inc/CaloHitMCNavigator.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"


#include "CLHEP/Vector/ThreeVector.h"


// C++ includes
#include <map>
#include <iostream>


namespace mu2e {



     class CaloClusterMCUtil {

	   public:

	       CaloClusterMCUtil(double edep, double time, double mom, CLHEP::Hep3Vector const& pos): 
	          _edepTot(edep), _time(time),_mom(mom),_pos(pos)
	          {};


	       void update(double edep, double time, double mom, CLHEP::Hep3Vector const& pos)
	       {
        	  _edepTot += edep;	
		  if ( mom > _mom ) {_time=time; _pos=pos;_mom=mom;}
	       }


	       double                   edepTot()  const {return _edepTot;}
	       double                   time()     const {return _edepTot;}
	       double                   momentum() const {return _mom;}
	       CLHEP::Hep3Vector const& position() const {return _pos;}


	   private:

	       double                _edepTot;
	       double                _time;
	       double                _mom;
	       CLHEP::Hep3Vector     _pos;
     };












    class CaloContentMC {


	 public:

             CaloContentMC(CaloHitMCNavigator const& caloHitMCNavi, CaloCluster const& cluster);
             ~CaloContentMC(){};
	     
	     std::vector<art::Ptr<SimParticle> > const& simBase()  const {return _simBase;}
	     std::vector<art::Ptr<SimParticle> > const& simPart()  const {return _simPart;}
	     std::vector<double>                 const& edepTot()  const {return _edep;}
	     std::vector<double>                 const& time()     const {return _time;}
	     std::vector<double>                 const& momentum() const {return _mom;}
	     std::vector<CLHEP::Hep3Vector>      const& position() const {return _pos;}
	     
             
	     bool hasConversion(); 



	 private:

             typedef art::Ptr<SimParticle> SimParticlePtr;
             typedef std::map<SimParticlePtr, CaloClusterMCUtil> SimUtilMap;

	     void fillCluster(); 
             void fillMaps(CaloHitSimPartMC const& hitSim, std::set<SimParticlePtr >& simBaseSet);


	     CaloHitMCNavigator const*      _navi;
	     CaloCluster const*             _cluster;
	     SimUtilMap                     _simUtilMap;
             std::vector<SimParticlePtr>    _simBase;

	     std::vector<SimParticlePtr>    _simPart;
	     std::vector<double>            _edep;
	     std::vector<double>            _time;
	     std::vector<double>            _mom;
	     std::vector<CLHEP::Hep3Vector> _pos;

    };


} 

#endif
