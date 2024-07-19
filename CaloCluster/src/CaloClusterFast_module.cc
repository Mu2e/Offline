#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"

#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"

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
            fhicl::Atom<art::InputTag>  caloHitCollection { Name("caloHitCollection"), Comment("Calo Hit collection")};
            fhicl::Atom<double>         EminSeed          { Name("EminSeed"),          Comment("Minimum energy for a hit to be a cluster seed") };
            fhicl::Atom<double>         EnoiseCut         { Name("EnoiseCut"),         Comment("Minimum energy for a hit to be in a cluster") };
            fhicl::Atom<double>         ExpandCut         { Name("ExpandCut"),         Comment("Minimum energy for a hit to expand cluster") };
            fhicl::Atom<double>         deltaTime         { Name("deltaTime"),         Comment("Maximum time difference between seed and hit in cluster") };
            fhicl::Atom<double>         timeOffset        { Name("timeOffset"),        Comment("Time offset to add to base cluster time") };
            fhicl::Atom<int>            minSiPMPerHit     { Name("minSiPMPerHit"),     Comment("Minimum number of SiPM contributing to the hit") };
            fhicl::Atom<bool>           extendSearch      { Name("extendSearch"),      Comment("Search next-next neighbors for clustering") };
            fhicl::Atom<int>            diagLevel         { Name("diagLevel"),         Comment("Diag level"),0 };
        };

        explicit CaloClusterFast(const art::EDProducer::Table<Config>& config) :
          EDProducer{config},
          caloHitToken_  {consumes<CaloHitCollection>(config().caloHitCollection())},
          EminSeed_      (config().EminSeed()),
          EnoiseCut_     (config().EnoiseCut()),
          ExpandCut_     (config().ExpandCut()),
          deltaTime_     (config().deltaTime()),
          timeOffset_     (config().timeOffset()),
          minSiPMPerHit_ (config().minSiPMPerHit()),
          extendSearch_  (config().extendSearch()),
          diagLevel_     (config().diagLevel())
        {
           produces<CaloClusterCollection>();
        }

        void produce(art::Event& e) override;


     private:
        art::ProductToken<CaloHitCollection> caloHitToken_;
        double            EminSeed_;
        double            EnoiseCut_;
        double            ExpandCut_;
        double            deltaTime_;
        double            timeOffset_;
        int               minSiPMPerHit_;
        bool              extendSearch_;
        int               diagLevel_;

        void makeClusters(CaloClusterCollection&, const art::Handle<CaloHitCollection>&);
        void fillCluster(const Calorimeter&, const art::Handle<CaloHitCollection>&, const CaloHitCollection&,
                         const std::vector<size_t>&, CaloClusterCollection&);
  };



  void CaloClusterFast::produce(art::Event& event)
  {
      art::Handle<CaloHitCollection> caloHitsHandle = event.getHandle<CaloHitCollection>(caloHitToken_);

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

      std::vector<size_t> hits;
      hits.reserve(caloHits.size());
      for (size_t i=0;i<caloHits.size();++i) if (caloHits[i].energyDep() > EnoiseCut_ && caloHits[i].nSiPMs() >= minSiPMPerHit_) hits.emplace_back(i);

      auto functorTime = [&caloHits,&hits](auto a, auto b) {return caloHits[a].time() < caloHits[b].time();};
      std::stable_sort(hits.begin(),hits.end(),functorTime);

      auto iterSeed = hits.begin();
      while (iterSeed != hits.end())
      {
          //find the first hit above the energy threshold, and the last hit within the required time window
          const CaloHit& hitSeed = caloHits[*iterSeed];
          if (*iterSeed==hits.size() || hitSeed.energyDep()< EminSeed_) {++iterSeed; continue;}
          double timeStart = hitSeed.time();

          //find the range around the seed time to search for other hits to form clusters
          auto iterStart(iterSeed), iterStop(iterSeed);
          while (iterStop  != hits.end()   && (*iterStop==hits.size()  || caloHits[*iterStop].time() - timeStart < deltaTime_))  ++iterStop;
          while (iterStart != hits.begin() && (*iterStart==hits.size() || timeStart - caloHits[*iterStart].time() < deltaTime_)) --iterStart;
          ++iterStart;

          //start the clustering algorithm for the hits between iStart and iStop
          std::queue<int> crystalToVisit;
          std::vector<bool> isVisited(cal.nCrystals());

          //put the first hit in the cluster list
          std::vector<size_t> clusterList{*iterSeed};
          auto seedId = caloHits[*iterSeed].crystalID();
          crystalToVisit.push(seedId);
          *iterSeed=hits.size();

          // loop until all seeds are processed
          while (!crystalToVisit.empty())
          {
              auto visitId = crystalToVisit.front();
              isVisited[visitId]=true;

              auto neighborsId = cal.crystal(visitId).neighbors();
              if (extendSearch_) std::copy(cal.nextNeighbors(visitId).begin(), cal.nextNeighbors(visitId).end(), std::back_inserter(neighborsId));

              for (const auto& iId : neighborsId)
              {
                  if (isVisited[iId]) continue;
                  isVisited[iId] = true;

                  //loop over the caloHits, check if one is in the neighbor list and add it to the cluster
                  for (auto it=iterStart; it != iterStop; ++it)
                  {
                      if (*it==hits.size()) continue;
                      if (caloHits[*it].crystalID() != iId) continue;

                      if (caloHits[*it].energyDep() > ExpandCut_) crystalToVisit.push(iId);
                      clusterList.push_back(*it);
                      *it = hits.size();
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
                                      const CaloHitCollection& caloHits, const std::vector<size_t>& clusterList,
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
      double time            = seedHit.time() + timeOffset_;
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

DEFINE_ART_MODULE(mu2e::CaloClusterFast)
