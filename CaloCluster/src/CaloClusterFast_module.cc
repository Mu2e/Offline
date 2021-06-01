#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"

#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "RecoDataProducts/inc/CaloHit.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"

#include <iostream>
#include <string>
#include <list>
#include <set>
#include <vector>
#include <queue>


namespace mu2e {

  class CaloClusterFast : public art::EDProducer
  {
     public:        
        struct Config
        {
            using Name    = fhicl::Name;
            using Comment = fhicl::Comment;
            fhicl::Atom<std::string>  caloCrystalModuleLabel{ Name("caloCrystalModuleLabel"), Comment("Calo Crystal module label")};
            fhicl::Atom<double>       EminSeed              { Name("EminSeed"),               Comment("Minimum energy for a hit to be a cluster seed") }; 
            fhicl::Atom<double>       EnoiseCut             { Name("EnoiseCut"),              Comment("Minimum energy for a hit to be in a cluster") }; 
            fhicl::Atom<double>       ExpandCut             { Name("ExpandCut"),              Comment("Minimum energy for a hit to expand cluster") }; 
            fhicl::Atom<double>       deltaTime             { Name("deltaTime"),              Comment("Maximum time difference between seed and hit in cluster") }; 
            fhicl::Atom<bool>         extendSearch          { Name("extendSearch"),           Comment("Search next-next neighbors for clustering") }; 
            fhicl::Atom<int>          diagLevel             { Name("diagLevel"),              Comment("Diag level"),0 }; 
        };

        explicit CaloClusterFast(const art::EDProducer::Table<Config>& config) :
          EDProducer{config},
          caloCrystalToken_{consumes<CaloHitCollection>(config().caloCrystalModuleLabel())},
          EminSeed_        (config().EminSeed()),
          EnoiseCut_       (config().EnoiseCut()),
          ExpandCut_       (config().ExpandCut()),
          deltaTime_       (config().deltaTime()),
          extendSearch_    (config().extendSearch()),
          diagLevel_       (config().diagLevel())
        {
           produces<CaloClusterCollection>();
        }

        void produce(art::Event& e) override;


     private:
        art::ProductToken<CaloHitCollection> caloCrystalToken_;
        double            EminSeed_;
        double            EnoiseCut_;
        double            ExpandCut_;
        double            deltaTime_;
        bool              extendSearch_;
        int               diagLevel_;

        void makeClusters(CaloClusterCollection&, const art::Handle<CaloHitCollection>&);
        void fillCluster(const Calorimeter&, const art::Handle<CaloHitCollection>&, const CaloHitCollection&,
                         const std::vector<int>&, CaloClusterCollection&);
  };



  void CaloClusterFast::produce(art::Event& event)
  {
      art::Handle<CaloHitCollection> caloHitsHandle = event.getHandle<CaloHitCollection>(caloCrystalToken_);
      
      auto caloClusters = std::make_unique<CaloClusterCollection>();
      makeClusters(*caloClusters,caloHitsHandle);

      event.put(std::move(caloClusters));
  }


  //----------------------------------------------------------------------------------------------------------
  void CaloClusterFast::makeClusters(CaloClusterCollection& caloClusters, const art::Handle<CaloHitCollection>& caloHitsHandle)
  {
      const Calorimeter& cal = *(GeomHandle<Calorimeter>());
      const CaloHitCollection& caloHits(*caloHitsHandle);
      if (caloHits.empty()) return;

      std::vector<int> hits;
      hits.reserve(caloHits.size());
      for (unsigned i=0;i<caloHits.size();++i) if (caloHits[i].energyDep() > EnoiseCut_) hits.emplace_back(i);      
      auto functorTime = [&caloHits](int a, int b) {return caloHits[a].time() < caloHits[b].time();};
      std::sort(hits.begin(),hits.end(),functorTime);

      auto iterSeed = hits.begin();
      while (iterSeed != hits.end())
      {
          //find the first hit above the energy threshold, and the last hit within the required time window
          const CaloHit& hitSeed = caloHits[*iterSeed];
          if (*iterSeed==-1 || hitSeed.energyDep()< EminSeed_) {++iterSeed; continue;}
          double timeStart = hitSeed.time();

          //find the range around the seed time to search for other hits to form clusters
          auto iterStart(iterSeed), iterStop(iterSeed);
          while (iterStop  != hits.end()   && (*iterStop==-1  || caloHits[*iterStop].time() - timeStart < deltaTime_))  ++iterStop;
          while (iterStart != hits.begin() && (*iterStart==-1 || timeStart - caloHits[*iterStart].time() < deltaTime_)) --iterStart;
          ++iterStart; 

          //start the clustering algorithm for the hits between iStart and iStop
          std::queue<int> crystalToVisit;
          std::vector<bool> isVisited(cal.nCrystal()); 

          //put the first hit in the cluster list
          std::vector<int> clusterList{*iterSeed};
          int seedId = caloHits[*iterSeed].crystalID();
          crystalToVisit.push(seedId);
          *iterSeed=-1;

          // loop until all seeds are processed
          while (!crystalToVisit.empty())
          {            
              int visitId = crystalToVisit.front();
              isVisited[visitId]=true;

              std::vector<int>  neighborsId = cal.crystal(visitId).neighbors();
              if (extendSearch_) std::copy(cal.nextNeighbors(visitId).begin(), cal.nextNeighbors(visitId).end(), std::back_inserter(neighborsId));

              for (const auto& iId : neighborsId)
              {
                  if (isVisited[iId]) continue;
                  isVisited[iId] = true;

                  //loop over the caloHits, check if one is in the neighbor list and add it to the cluster
                  for (auto it=iterStart; it != iterStop; ++it)
                  {
                      if (*it==-1) continue;
                      if (caloHits[*it].crystalID() != iId) continue;

                      if (caloHits[*it].energyDep() > ExpandCut_) crystalToVisit.push(iId);
                      clusterList.push_back(*it);
                      *it = -1;
                  }
               }
               crystalToVisit.pop();                 
           }

           auto functorEnergy = [&caloHits](int a, int b) {return caloHits[a].energyDep() > caloHits[b].energyDep();};
           std::sort(clusterList.begin(),clusterList.end(),functorEnergy);

           fillCluster(cal, caloHitsHandle, caloHits, clusterList, caloClusters );

           ++iterSeed;
        }      
   }

   
  //----------------------------------------------------------------------------------------------------------
  void CaloClusterFast::fillCluster(const Calorimeter& cal, const art::Handle<CaloHitCollection>& caloHitsHandle, 
                                      const CaloHitCollection& caloHits, const std::vector<int>& clusterList, 
                                      CaloClusterCollection& caloClusters)
  {
      std::vector<art::Ptr<CaloHit>> caloHitsPtrs;
      double totalEnergy(0), xCOG(0), yCOG(0);
 
      for (auto idx : clusterList) 
      {
          int crID        = caloHits[idx].crystalID();
          totalEnergy    += caloHits[idx].energyDep();
          xCOG           += cal.crystal(crID).localPosition().x()*caloHits[idx].energyDep();
          yCOG           += cal.crystal(crID).localPosition().y()*caloHits[idx].energyDep(); 
          caloHitsPtrs.push_back(art::Ptr<CaloHit>(caloHitsHandle,idx));
      }

      xCOG /= totalEnergy;
      yCOG /= totalEnergy;

      const CaloHit& seedHit = caloHits[clusterList[0]];
      double time            = seedHit.time();
      int    iDisk           = cal.crystal(seedHit.crystalID()).diskID();
 
      caloClusters.emplace_back(CaloCluster(iDisk,time,0.0,totalEnergy,0.0,CLHEP::Hep3Vector(xCOG,yCOG,0),
                                            caloHitsPtrs,clusterList.size(),false));

      if (diagLevel_ > 0)
      {
          std::cout<<"This cluster contains "<<clusterList.size()<<" crystals, id= ";
          for (auto val : clusterList) std::cout<<caloHits[val].crystalID()<<" ";
          std::cout<<" with energy="<<totalEnergy<<" and time="<<time<<std::endl;;
      }
  }

}

DEFINE_ART_MODULE(mu2e::CaloClusterFast);
