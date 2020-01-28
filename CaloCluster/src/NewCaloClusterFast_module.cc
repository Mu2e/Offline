//Author: S Middleton
//Date: Nov 2019
//Purpose: For the purpose of fast clustering from crystal hits


#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"

#include "GeometryService/inc/GeomHandle.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/CalorimeterCalibrations.hh"
#include "RecoDataProducts/inc/CaloDigi.hh"
#include "RecoDataProducts/inc/CaloDigiCollection.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHit.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"

#include <iostream>
#include <string>
#include <queue>

namespace mu2e {
	class NewCaloClusterFast : public art::EDProducer {

		public:

			explicit NewCaloClusterFast(fhicl::ParameterSet const& pset) :
			art::EDProducer{pset},
			caloCrystalToken_{consumes<CaloCrystalHitCollection>(pset.get<std::string>("caloCrystalModuleLabel"))},
			extendSecond_(pset.get<bool>("extendSecond")),
			minEnergy_( pset.get<double>("minEnergy")),
			diagLevel_(pset.get<int>("diagLevel",0)),
			deltaTime_(pset.get<double>("deltaTime"))
			{
				produces<CaloClusterCollection>();
			}
			typedef std::list<const CaloCrystalHit*>    CaloCrystalList;
			virtual ~NewCaloClusterFast() {};
			virtual void produce(art::Event& e) override;

		private:
			art::ProductToken<CaloCrystalHitCollection> const caloCrystalToken_;

			bool         extendSecond_;
			double       minEnergy_;
			int          diagLevel_;
			double       deltaTime_;
			art::ProductID                _crystalHitsPtrID;
			art::EDProductGetter const*  _crystalHitsPtrGetter;

			void MakeClusters(CaloClusterCollection& recoClusters, const art::Handle<CaloCrystalHitCollection> & recoCrystalHit);
  
	};

	void NewCaloClusterFast::produce(art::Event& event)
 	{
		if (diagLevel_ > 0) std::cout<<"[NewCaloClusterFast::produce] begin"<<std::endl;

		art::Handle<CaloCrystalHitCollection> caloCrystalHitsHandle;
		bool const success = event.getByToken(caloCrystalToken_, caloCrystalHitsHandle);
		if (!success) return;

		auto recoClustersColl = std::make_unique<CaloClusterCollection>();
		recoClustersColl->reserve(10);

		_crystalHitsPtrID     = event.getProductID<CaloCrystalHitCollection>();
		_crystalHitsPtrGetter = event.productGetter(_crystalHitsPtrID);

		MakeClusters( *recoClustersColl, caloCrystalHitsHandle);

		if ( diagLevel_ > 3 )
		{
		   printf("[NewCaloClusterFast::produce] produced RecoCrystalHits ");
		   printf(", recoClustersColl size  = %i \n", int(recoClustersColl->size()));
		}
		event.put(std::move(recoClustersColl));
		if (diagLevel_ > 0) std::cout<<"[NewCaloClusterFast::produce] end"<<std::endl;
		return;
	}

  	void NewCaloClusterFast::MakeClusters(CaloClusterCollection& recoClusters, const art::Handle<CaloCrystalHitCollection> & CaloCrystalHitsHandle)
  	{
		//Find the hit collection 
		const CaloCrystalHitCollection& recoCrystalHits(*CaloCrystalHitsHandle);

		//Find calorimter:
		mu2e::GeomHandle<mu2e::Calorimeter> ch;
		const Calorimeter* cal = ch.get();

		if (recoCrystalHits.empty()) return;

		std::vector<CaloCrystalList>      clusterList, caloIdHitMap(cal->nCrystal());
		std::list<const CaloCrystalHit*>  seeds_;

		//Loop through collection, fill seeds and map:
		for (const auto& hit : recoCrystalHits)
		{
			if (hit.energyDep() < minEnergy_) continue;
			caloIdHitMap[hit.id()].push_back(&hit);
			seeds_.push_back(&hit);
		}

		//Sort seeds interms of energy deposited.
		seeds_.sort([](const CaloCrystalHit* a, const CaloCrystalHit* b) {return a->energyDep() > b->energyDep();});
		//check the seeds arnt empty, if not proceed:
		while( !seeds_.empty() )
		{
			//find seed:
			const CaloCrystalHit* crystalSeed = *seeds_.begin();
			//check it fulfils criteria on energy:
			if (crystalSeed->energyDep() < minEnergy_) break;
			//make some lists for cluster finding:
			CaloCrystalList List;
			std::queue<int> crystalToVisit_;

			double seedTime_ = crystalSeed->time();
			//Add seed to front of this List of crystals:
			List.push_front(crystalSeed); 
			//Add crystal ID to list of crytals to visit:

			crystalToVisit_.push(crystalSeed->id());  
			if (extendSecond_) for (const auto& nneighbor : cal->nextNeighbors(crystalSeed->id())) crystalToVisit_.push(nneighbor);
			//Find the crystal from map, which is assoiciated with this seed:
			CaloCrystalList& cryL = caloIdHitMap[crystalSeed->id()];
			cryL.erase(std::find(cryL.begin(), cryL.end(), crystalSeed));
			std::vector<bool>      isVisited_;
			isVisited_.reserve(crystalToVisit_.size());

			while (!crystalToVisit_.empty())
            		{            
 		        
				 int visitId         = crystalToVisit_.front();
				 isVisited_[visitId] = 1;
				 //loop through the neighbours if the crystal and check:
				 std::vector<int> const& neighborsId = cal->crystal(visitId).neighbors();
				 for (auto& iId : neighborsId)
				 {               
					//if visited --> skip
					if (isVisited_[iId]) continue;
					//if not then assign a flag to say its been done now:
					isVisited_[iId]=1;
					//copy list:
					CaloCrystalList& list = caloIdHitMap[iId];
					//iterate:
					auto it=list.begin();
					while(it != list.end())
                     			{
						//get hit:
						 CaloCrystalHit const* hit = *it;
							 //time check:
						 if (std::abs(hit->time() - seedTime_) < deltaTime_) 
						{ 
						     	//check energy:
							if (hit->energyDep() < minEnergy_) { crystalToVisit_.push(iId); }
							 //add hit:
							List.push_front(hit);
							//remove:
							it = list.erase(it);   
						} 
						else {++it;}
             				} 
                     
                 		}	
                                       
                 		crystalToVisit_.pop();                 
            		} //ends while 
			//sort the List of crystals again interms of energy:
			List.sort([] (CaloCrystalHit const* lhs, CaloCrystalHit const* rhs) {return lhs->energyDep() > rhs->energyDep();} );               
			//add list to lis-of-lists for cluster making:
			clusterList.push_back(List);
			//remove from seeds:
			for (const auto& hit : List) { seeds_.remove(hit); } 
     		}//ends seeds
      
		//loop through list for clustering:
		// Note: ClusterList is a list of crystal_lists, each element of the cluster list is a cluster
		for (auto cluster : clusterList) { 
        
			const CaloCrystalHit* caloCrystalHitBase = &recoCrystalHits.front();
			    //create pointer --> this will be associated with the cluster if include set
			std::vector<art::Ptr<CaloCrystalHit>> caloCrystalHitsPtrVector;
			    //need to add up these over all crystals in cluster:
			double totalEnergy(0),totalEnergyErr(0), xcl(0), ycl(0), ncry(0);
			    //loop:
			for (auto clusterPrt : cluster) //Note: a clusterPrt is a crystal hit
			{
				int    crId = clusterPrt->id();
				totalEnergy    += clusterPrt->energyDep();
				totalEnergyErr += clusterPrt->energyDepErr()*clusterPrt->energyDepErr();
				xcl += cal->crystal(crId).localPosition().x()*clusterPrt->energyDep();
				ycl += cal->crystal(crId).localPosition().y()*clusterPrt->energyDep();
				size_t idx = (clusterPrt - caloCrystalHitBase);
				caloCrystalHitsPtrVector.push_back( art::Ptr<CaloCrystalHit>(CaloCrystalHitsHandle,idx) ); 
				ncry++;
			}
			//make cluster:
			totalEnergyErr = sqrt(totalEnergyErr);
			xcl = xcl/totalEnergy;
			ycl = ycl/totalEnergy;

			double time    = (*cluster.begin())->time();
			double timeErr = (*cluster.begin())->timeErr();
			const auto& seed  = **caloCrystalHitsPtrVector.begin();
			int iSection = cal->crystal(seed.id()).diskId(); 

			CaloCluster Endcluster(iSection,time,timeErr,totalEnergy,totalEnergyErr,caloCrystalHitsPtrVector,ncry,0.0);
			//Get centre-of-g:
			Endcluster.cog3Vector(CLHEP::Hep3Vector(xcl,ycl,0));
			//Add cluster to event:
			recoClusters.emplace_back(std::move(Endcluster));

		} //ends clusters
 	}//end function
}
DEFINE_ART_MODULE(mu2e::NewCaloClusterFast);
