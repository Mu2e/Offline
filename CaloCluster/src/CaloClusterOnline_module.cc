// Author: S Middleton
// Date: Feb 2020
// Purpose: Fast Online clustering
// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "cetlib_except/exception.h"

#include "CaloCluster/inc/ClusterFinder.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "RecoDataProducts/inc/CaloHit.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"

// Other includes.
#include <iostream>
#include <string>
#include <list>
#include <vector>


namespace mu2e {
    class CaloClusterOnline : public art::EDProducer {
        
        public:

            typedef std::vector<const CaloHit*>  CaloCrystalVec;
            typedef std::list<const CaloHit*>    CaloCrystalList;


            explicit CaloClusterOnline(fhicl::ParameterSet const& pset) :
            art::EDProducer{pset},
            caloCrystalToken_{consumes<CaloHitCollection>(pset.get<std::string>("caloCrystalModuleLabel"))},
            minClusterEnergy_(pset.get<double>("minClusterEnergy",0)),
            EminSeed_(pset.get<double>("EminSeed")),
            EnoiseCut_(pset.get<double>("EnoiseCut")),
            ExpandCut_(pset.get<double>("ExpandCut")),
            timeCut_(pset.get<double>("timeCut")),
            deltaTime_(pset.get<double>("deltaTime")),
            diagLevel_(pset.get<int>("diagLevel",0)),
            messageCategory_("CLUSTER")
            {
                produces<CaloClusterCollection>();
            }

            void produce(art::Event& e) override;

        private:

            art::ProductToken<CaloHitCollection> const caloCrystalToken_;
            double            minClusterEnergy_;
            double            EminSeed_;
            double            EnoiseCut_;
            double            ExpandCut_;
            double            timeCut_;
            double            deltaTime_;
            int               diagLevel_;
            const std::string messageCategory_;

            void MakeOnlineClusters(CaloClusterCollection& caloClusters,
                       const art::Handle<CaloHitCollection>& CaloHitsHandle);

            void FillOnlineCluster(CaloClusterCollection& caloClustersColl, const CaloCrystalList& clusterList,
                 const art::Handle<CaloHitCollection>& CaloHitsHandle, const Calorimeter& cal);
    };


    void CaloClusterOnline::produce(art::Event& event)
    {
        if (diagLevel_ > 0) std::cout<<"[CaloClusterOnlines::produce] begin"<<std::endl;

        // Check that calorimeter geometry description exists
        art::ServiceHandle<GeometryService> geom;
        if( !(geom->hasElement<Calorimeter>()) ) return;

        // Get handles to calorimeter crystal hits
        art::Handle<CaloHitCollection> CaloHitsHandle;
        bool const success = event.getByToken(caloCrystalToken_, CaloHitsHandle);
        if (!success) return;
        auto recoClustersColl = std::make_unique<CaloClusterCollection>();
        if ( diagLevel_ > 0 )
        {
            std::cout<<"[OnlineClusterMaker::produce] No. RecoCrystalHits: "<<recoClustersColl->size()<<std::endl;
        }
        MakeOnlineClusters(*recoClustersColl,CaloHitsHandle);
        event.put(std::move(recoClustersColl));

    if (diagLevel_ > 0) std::cout<<"[CaloClusterOnline::produce] end"<<std::endl;
    return;
  }
    
  void CaloClusterOnline::MakeOnlineClusters(CaloClusterCollection& recoClusters, const art::Handle<CaloHitCollection> & CaloHitsHandle)
      {
        const Calorimeter& cal = *(GeomHandle<Calorimeter>());
        const CaloHitCollection& CaloHits(*CaloHitsHandle);
        if (CaloHits.empty()) return;

        std::vector<CaloCrystalList>      clusterList, caloIdHitMap(cal.nCrystal());
        std::list<const CaloHit*>  seedList;
        
        for (const auto& hit : CaloHits)
        {
            if (hit.energyDep() <  EnoiseCut_) continue;
            caloIdHitMap[hit.id()].push_back(&hit);
            seedList.push_back(&hit);
        }

        seedList.sort([](const CaloHit* a, const CaloHit* b) {return a->energyDep() > b->energyDep();});

        while( !seedList.empty() )
        {
            const CaloHit* crystalSeed = *seedList.begin();
            if (crystalSeed->energyDep() < EminSeed_) break;

            ClusterFinder finder(cal,crystalSeed,deltaTime_, ExpandCut_, true);
            finder.formCluster(caloIdHitMap);
	    
	    clusterList.push_back(finder.clusterList());
            for (const auto& hit: finder.clusterList()) seedList.remove(hit);

      }

        for (auto cluster : clusterList)  FillOnlineCluster(recoClusters,
						cluster,CaloHitsHandle,cal);
  }

  void CaloClusterOnline::FillOnlineCluster(CaloClusterCollection& caloClustersColl, const CaloCrystalList& clusterPtrList, const art::Handle<CaloHitCollection>& CaloHitsHandle, const Calorimeter& cal)
  {

    const CaloHitCollection& recoCrystalHits(*CaloHitsHandle);
    const CaloHit* CaloHitBase = &recoCrystalHits.front();
    std::vector<art::Ptr<CaloHit>> caloHitsPtrVector;

    double totalEnergy(0),totalEnergyErr(0), xcl(0), ycl(0), ncry(0);

    for (auto clusterPrt : clusterPtrList)
    {
        int    crId = clusterPrt->id();
        totalEnergy    += clusterPrt->energyDep();
        totalEnergyErr += clusterPrt->energyDepErr()*clusterPrt->energyDepErr();
        xcl += cal.crystal(crId).localPosition().x()*clusterPrt->energyDep();
        ycl += cal.crystal(crId).localPosition().y()*clusterPrt->energyDep();
        size_t idx = (clusterPrt - CaloHitBase);
        caloHitsPtrVector.push_back( art::Ptr<CaloHit>(CaloHitsHandle,idx) );
        ncry++;
       
    }

    if (totalEnergy < minClusterEnergy_)  return;

    totalEnergyErr = sqrt(totalEnergyErr);
    xcl = xcl/totalEnergy;
    ycl = ycl/totalEnergy;

    double time    = (*clusterPtrList.begin())->time();
    double timeErr = (*clusterPtrList.begin())->timeErr();
    const auto& seed  = **caloHitsPtrVector.begin();
    int iSection = cal.crystal(seed.id()).diskId();

    CaloCluster Endcluster(iSection,time,timeErr,totalEnergy,
				totalEnergyErr,caloHitsPtrVector,ncry,0.0);
    Endcluster.cog3Vector(CLHEP::Hep3Vector(xcl,ycl,0));
    caloClustersColl.emplace_back(Endcluster);
  }
}

DEFINE_ART_MODULE(mu2e::CaloClusterOnline);


