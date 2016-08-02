//
// Module to produce the calorimeter clusters from proto-clusters. See MakeCaloProtocluster for the proto-cluster formation.
//  
// The strategy is to attach the split-off to the main cluster first, take the closest cluster if the split-off time is compatible 
// with several clusters. Then associate energetic clusters between them, including the reattached split-off of each cluster in the 
// comparison. 
//
// Note: the cluster center-of-gravity is calculated in the calorimeter section front face frame
//
// Original author: B. Echenard


// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


// Mu2e includes.
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CaloCluster/inc/ClusterAssociator.hh"
#include "CaloCluster/inc/ClusterMoments.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "RecoDataProducts/inc/CaloProtoClusterCollection.hh"

// Other includes.
#include "cetlib/exception.h"

#include "TH1D.h"
#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <memory>
#include <tuple>
#include <array>



namespace mu2e {



class CaloClusterFromProtoCluster : public art::EDProducer {


     public:

             explicit CaloClusterFromProtoCluster(fhicl::ParameterSet const& pset) :
             _diagLevel(pset.get<int>("diagLevel",0)),
             _caloClusterModuleLabel(pset.get<std::string>("caloClusterModuleLabel")),
             _producerNameMain(pset.get<std::string>("mainClusterCollName")),
             _producerNameSplit(pset.get<std::string>("splitClusterCollName")),
             _deltaTimePlus(pset.get<double>("deltaTimePlus")),    
             _deltaTimeMinus(pset.get<double>("deltaTimeMinus")),  
             _maxDistSplit(pset.get<double>("maxDistSplit")),                
             _maxDistMain(pset.get<double>("maxDistMain")),                
             _cogTypeName(pset.get<std::string>("cogTypeName")),
             _cogType(ClusterMoments::Linear),                
             _messageCategory("CLUSTER")
             {
                 produces<CaloClusterCollection>();
             }

             virtual ~CaloClusterFromProtoCluster() { }
             virtual void beginJob();
             void produce( art::Event& e);


     private:

             int                _diagLevel;
             std::string        _caloClusterModuleLabel;
             std::string        _producerNameMain;
             std::string        _producerNameSplit;
             double             _deltaTimePlus;
             double             _deltaTimeMinus;
             double             _maxDistSplit;
             double             _maxDistMain;
             std::string        _cogTypeName;             
             ClusterMoments::cogtype _cogType;
             const std::string  _messageCategory;
            
             void makeCaloClusters(CaloClusterCollection& caloClusters, 
                                   CaloProtoClusterCollection const& caloClustersMain, 
                                   CaloProtoClusterCollection const& caloClustersSplit);
     
             std::tuple<double,double,double> CalcEnergyLayer(Calorimeter const & cal, 
                                              std::vector<art::Ptr<CaloCrystalHit>> const& caloCrystalHitsPtrVector);
     
     };



     void CaloClusterFromProtoCluster::beginJob()
     {
         if (_cogTypeName.compare("Linear"))    _cogType = ClusterMoments::Linear;         
         if (_cogTypeName.compare("Logarithm")) _cogType = ClusterMoments::Logarithm;         
     }




     void CaloClusterFromProtoCluster::produce(art::Event& event) {


         // Check that calorimeter geometry description exists
         art::ServiceHandle<GeometryService> geom;
         if( !(geom->hasElement<Calorimeter>()) ) return;

         //Get calo cluster
         art::Handle<CaloProtoClusterCollection> caloClustersMainHandle;
         event.getByLabel(_caloClusterModuleLabel, _producerNameMain, caloClustersMainHandle);
         CaloProtoClusterCollection const& caloClustersMain(*caloClustersMainHandle);

         art::Handle<CaloProtoClusterCollection> caloClustersSplitHandle;
         event.getByLabel(_caloClusterModuleLabel, _producerNameSplit, caloClustersSplitHandle);
         CaloProtoClusterCollection const& caloClustersSplit(*caloClustersSplitHandle);

         //Create a new CaloCluster collection and fill it
         std::unique_ptr<CaloClusterCollection> caloClusters(new CaloClusterCollection);
         makeCaloClusters(*caloClusters, caloClustersMain, caloClustersSplit);

         event.put(std::move(caloClusters));

     }





     void CaloClusterFromProtoCluster::makeCaloClusters(CaloClusterCollection& caloClusters, 
                                       CaloProtoClusterCollection const& caloClustersMain, 
                                       CaloProtoClusterCollection const& caloClustersSplit)
     {

	Calorimeter const & cal = *(GeomHandle<Calorimeter>());

	// intermediate buffer for storing main clusters + split offs
	CaloProtoClusterCollection caloProtoClustersTemp;

	ClusterAssociator associator(cal);
	associator.associateSplitOff(caloClustersMain, caloClustersSplit, _deltaTimePlus, _deltaTimeMinus,_maxDistSplit);


	// associate split-off to main cluster
	for (unsigned int imain=0;imain<caloClustersMain.size();++imain)
	{ 
            std::vector<art::Ptr<CaloCrystalHit>> caloCrystalHitsPtrVector;
            bool isSplit(false);
            double totalEnergy(0);

            auto const& main = caloClustersMain.at(imain).caloCrystalHitPtrVector();
            auto const& seed = **main.begin();
            double seed_time = seed.time();


            //look at the main cluster
            if (_diagLevel > 1) std::cout<<"Associated main cluster "<<imain <<" at time "<<seed_time<<" with id= ";

            for (auto il = main.begin(); il !=main.end(); ++il)
            {
        	totalEnergy += (*il)->energyDep();
        	caloCrystalHitsPtrVector.push_back(*il);
        	if (_diagLevel > 1 ) std::cout<<(*il)->id()<<" "; 
            }


            //adding split-off here                  
            for (unsigned int isplit=0;isplit<caloClustersSplit.size();++isplit)
            {                    
        	if (associator.associatedSplitId(isplit) != imain) continue;
        	isSplit = true;

        	if (_diagLevel > 1) std::cout<<"with split-off with id= ";
        	auto const& split = caloClustersSplit.at(isplit).caloCrystalHitPtrVector();
        	for (auto il = split.begin(); il !=split.end(); ++il)
        	{
                    totalEnergy += (*il)->energyDep();
                    caloCrystalHitsPtrVector.push_back(*il);
                    if (_diagLevel > 2 ) std::cout<<(*il)->id()<<" "; 
        	}
            }
            if (_diagLevel > 1) std::cout<<std::endl;

            CaloProtoCluster CaloProtoCluster(seed_time,totalEnergy,caloCrystalHitsPtrVector,isSplit);
            caloProtoClustersTemp.push_back(CaloProtoCluster);
	}        




	// combine main clusters together (split-off included in main clusters at this point)
	associator.associateMain(caloProtoClustersTemp, _deltaTimePlus, _maxDistMain);


	//finally, form final clusters
	std::vector<int> flagProto(caloClustersMain.size(),0);
	for (unsigned int iproto=0;iproto<caloProtoClustersTemp.size();++iproto)
	{ 
            if (flagProto[iproto]) continue;                                    

            auto   caloCrystalHitsPtrVector = caloProtoClustersTemp.at(iproto).caloCrystalHitPtrVector();
            bool   isSplit                    = caloProtoClustersTemp.at(iproto).isSplit();
            double totalEnergy                = caloProtoClustersTemp.at(iproto).energyDep();

            for (int iassoc : associator.associatedMainId(iproto))
            {
        	if (_diagLevel > 1) std::cout<<"Associated to main cluster id="<<iproto<<"   main split="<<iassoc<<std::endl;
        	flagProto[iassoc] = 1;
        	totalEnergy += caloProtoClustersTemp.at(iassoc).energyDep();
        	isSplit = true;
        	caloCrystalHitsPtrVector.insert(caloCrystalHitsPtrVector.end(), 
                                               caloProtoClustersTemp.at(iassoc).caloCrystalHitPtrVector().begin(), 
                                               caloProtoClustersTemp.at(iassoc).caloCrystalHitPtrVector().end());
            }


            std::sort(caloCrystalHitsPtrVector.begin(),caloCrystalHitsPtrVector.end(),
                      [] (art::Ptr<CaloCrystalHit> const& lhs, art::Ptr<CaloCrystalHit> const& rhs) 
                          {return lhs->energyDep() > rhs->energyDep();} );


            auto const& seed   = **caloCrystalHitsPtrVector.begin();
            double seed_time   = seed.time();
            int seed_section   = cal.crystal(seed.id()).sectionId();


            if (_diagLevel > 2)
            {
        	std::cout<<"Making a new cluster with id= ";
        	for (auto il = caloCrystalHitsPtrVector.begin(); il !=caloCrystalHitsPtrVector.end(); ++il) std::cout<<(*il)->id()<<" "; 
            }


            CaloCluster caloCluster(seed_section,seed_time,totalEnergy,caloCrystalHitsPtrVector,isSplit);              

            ClusterMoments cogCalculator(cal,caloCluster, seed_section);
            cogCalculator.calculate(_cogType);
            caloCluster.cog3Vector(cogCalculator.cog());
            caloCluster.secondMoment(cogCalculator.secondMoment());
            caloCluster.angle(cogCalculator.angle());

             auto EnerLayer = CalcEnergyLayer(cal,caloCrystalHitsPtrVector);                
            caloCluster.energyRing(std::get<0>(EnerLayer),std::get<1>(EnerLayer),std::get<2>(EnerLayer));

            caloClusters.push_back(caloCluster);
	}
     }




     //----------------------------------------------------------------------------------------------------
     std::tuple<double,double,double> CaloClusterFromProtoCluster::CalcEnergyLayer(Calorimeter const & cal, 
                                      std::vector<art::Ptr<CaloCrystalHit>> const& caloCrystalHitsPtrVector) 
     {
         int seedId              = caloCrystalHitsPtrVector[0]->id();
         double seedEnergy       = caloCrystalHitsPtrVector[0]->energyDep();
         auto const neighborsId  = cal.crystal(seedId).neighbors();
         auto const nneighborsId = cal.crystal(seedId).nextNeighbors();
         
         double e1(seedEnergy),e9(seedEnergy),e25(seedEnergy);
         for (auto const& il : caloCrystalHitsPtrVector)
         {
            int crid = il->id();
            for (auto const& it : neighborsId)  if (it==crid) {e9 += il->energyDep(); e25 += il->energyDep(); break;}
            for (auto const& it : nneighborsId) if (it==crid) {e25 += il->energyDep(); break;}
         }
         std::tuple<double,double,double> res(e1,e9,e25);
         return res;
     }





}



using mu2e::CaloClusterFromProtoCluster;
DEFINE_ART_MODULE(CaloClusterFromProtoCluster);




