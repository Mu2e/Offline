//
// Original author B.Echenard
//

#ifndef CaloCluster_CaloContentMC_HH_
#define CaloCluster_CaloContentMC_HH_


// Mu2e includes
#include "CaloCluster/inc/CaloContentMC.hh"
#include "HitMakers/inc/CaloReadoutUtilities.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/CaloHitSimPartMCCollection.hh"
#include "Mu2eUtilities/inc/CaloHitMCNavigator.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "HitMakers/inc/CaloHitSimUtil.hh"


// C++ includes
#include <map>
#include <iostream>


namespace mu2e {


    class CaloContentMC{


	 public:

             CaloContentMC(CaloHitMCNavigator const& caloHitMCNavi, CaloCluster const& cluster);
             ~CaloContentMC(){};
	     
	     std::vector<art::Ptr<SimParticle> > const& simBase() const    {return _simBase;}
	     std::vector<art::Ptr<SimParticle> > const& simPart() const    {return _simPart;}
	     std::vector<art::Ptr<StepPointMC> > const& stepPoints() const {return _steps;}
	     std::vector<double>                 const& edepTot() const    {return _edep;}
             
	     bool hasConversion(); 



	 private:

             typedef art::Ptr<SimParticle> SimParticlePtr;
             typedef art::Ptr<StepPointMC> StepPointMCPtr;
             typedef std::map<SimParticlePtr, CaloHitSimUtil> SimStepMap;

	     CaloHitMCNavigator const*    _navi;
	     CaloCluster const*           _cluster;

	     SimStepMap _simStep;
	     std::set<SimParticlePtr >    _simBaseSet;
	     std::vector<SimParticlePtr > _simPart;
             std::vector<SimParticlePtr > _simBase;
	     std::vector<StepPointMCPtr > _steps;	     
	     std::vector<double>          _edep;

	     void fillCluster(); 
             void fillMaps(CaloHitSimPartMC const& hitSim);

    };


} // end namespace mu2e

#endif
