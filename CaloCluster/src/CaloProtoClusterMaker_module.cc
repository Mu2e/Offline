//
// This module produces simply connected clusters (aka proto-clusters) in two steps
//
// 1. Form Energetic proto-clusters
//    - start from a seed with an energy greater than some threshold and ad all simply connected hits
//      (simply connected = any two hits in a cluster can be joined by a continuous path of clusters in the crystal)
//    - repeat until all "energetic" seeds are exhausted
//
// 2. Form Split-offs: some clusters might have low-energy split-offs, and we need to find them
//    - filter the remaining unassigned hits to retain only those compatible with the time of the main clusters
//    - update the seed lists and form all remaining clusters as before
//
// Note 1: Seed do not need to be ordered by energy
// Note 2: The cluster time is taken as that of the most energetic hit -> potential for improvement (have fun)
// Note 3: Several optimization obscured the code for little gain, so I sticked to simplicity
//

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"

#include "Offline/CaloCluster/inc/ClusterFinder.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/CaloProtoCluster.hh"

#include <iostream>
#include <string>
#include <list>
#include <set>
#include <vector>


namespace mu2e {

  class CaloProtoClusterMaker : public art::EDProducer
  {
     public:
        typedef std::vector<const CaloHit*>  CaloCrystalVec;
        typedef std::list<const CaloHit*>    CaloCrystalList;

        struct Config
        {
            using Name    = fhicl::Name;
            using Comment = fhicl::Comment;
            fhicl::Atom<art::InputTag>  caloHitCollection  { Name("caloHitCollection"), Comment("CaloHit collection name")};
            fhicl::Atom<double>         EminSeed           { Name("EminSeed"),          Comment("Minimum energy for a hit to be a cluster seed") };
            fhicl::Atom<double>         EnoiseCut          { Name("EnoiseCut"),         Comment("Minimum energy for a hit to be in a cluster") };
            fhicl::Atom<double>         ExpandCut          { Name("ExpandCut"),         Comment("Minimum energy for a hit to expand cluster") };
            fhicl::Atom<bool>           addSecondRing      { Name("addSecondRing"),     Comment("Add secondary ring around crystal when forming clusters") };
            fhicl::Atom<double>         deltaTime          { Name("deltaTime"),         Comment("Maximum time difference between seed and hit in cluster") };
            fhicl::Atom<int>            diagLevel          { Name("diagLevel"),         Comment("Diag level"),0 };
        };

        explicit CaloProtoClusterMaker(const art::EDProducer::Table<Config>& config) :
          EDProducer{config},
          caloCrystalToken_{consumes<CaloHitCollection>(config().caloHitCollection())},
          EminSeed_        (config().EminSeed()),
          EnoiseCut_       (config().EnoiseCut()),
          ExpandCut_       (config().ExpandCut()),
          addSecondRing_   (config().addSecondRing()),
          deltaTime_       (config().deltaTime()),
          diagLevel_       (config().diagLevel())
        {
           produces<CaloProtoClusterCollection>("main");
           produces<CaloProtoClusterCollection>("split");
        }

        void produce(art::Event& e) override;

     private:
        art::ProductToken<CaloHitCollection> caloCrystalToken_;
        double                               EminSeed_;
        double                               EnoiseCut_;
        double                               ExpandCut_;
        bool                                 addSecondRing_;
        double                               deltaTime_;
        int                                  diagLevel_;

        void makeProtoClusters (CaloProtoClusterCollection&,CaloProtoClusterCollection&, const art::Handle<CaloHitCollection>&);
        void filterByTime      (CaloCrystalList&, const std::vector<double>&);
        void fillCluster       (CaloProtoClusterCollection&, const CaloCrystalList&,const art::Handle<CaloHitCollection>&);
        void dump              (const std::string&, const std::vector<CaloCrystalList>&, std::set<const CaloHit*>);
  };


  void CaloProtoClusterMaker::produce(art::Event& event)
  {
      art::Handle<CaloHitCollection> CaloHitsHandle = event.getHandle<CaloHitCollection>(caloCrystalToken_);

      auto caloProtoClustersMain  = std::make_unique<CaloProtoClusterCollection>();
      auto caloProtoClustersSplit = std::make_unique<CaloProtoClusterCollection>();
      makeProtoClusters(*caloProtoClustersMain,*caloProtoClustersSplit,CaloHitsHandle);

      event.put(std::move(caloProtoClustersMain),  "main");
      event.put(std::move(caloProtoClustersSplit), "split");
  }


  //----------------------------------------------------------------------------------------------------------
  void CaloProtoClusterMaker::makeProtoClusters(CaloProtoClusterCollection& caloProtoClustersMain,
                                                CaloProtoClusterCollection& caloProtoClustersSplit,
                                                const art::Handle<CaloHitCollection> & CaloHitsHandle)
  {
      const Calorimeter& cal = *(GeomHandle<Calorimeter>());
      const CaloHitCollection& CaloHits(*CaloHitsHandle);
      if (CaloHits.empty()) return;


      //declare and fill the hash map crystal_id -> list of CaloHits
      std::vector<CaloCrystalList> mainClusterList, splitClusterList, caloIdHitMap(cal.nCrystals());
      std::set<const CaloHit*> seedList;
      std::vector<double>      clusterTime;

      //fill data structures
      for (const auto& hit : CaloHits)
      {
          if (hit.energyDep() < EnoiseCut_) continue;
          caloIdHitMap[hit.crystalID()].push_back(&hit);
          if (hit.energyDep() > EminSeed_ ) seedList.insert(&hit);
      }

      if (diagLevel_ > 2) dump("Init", caloIdHitMap,seedList);



      //produce main clusters
      while( !seedList.empty() )
      {
          const CaloHit* crystalSeed = *seedList.begin();
          ClusterFinder finder(cal, crystalSeed, deltaTime_, ExpandCut_, addSecondRing_);
          finder.formCluster(caloIdHitMap);

          mainClusterList.push_back(finder.clusterList());
          clusterTime.push_back(crystalSeed->time());

          for (const auto& hit: finder.clusterList()) seedList.erase(hit);
      }


      //filter unneeded hits and fill new seeds
      for (unsigned i=0; i < caloIdHitMap.size(); ++i ) filterByTime(caloIdHitMap[i], clusterTime);
      for (const auto& liste: caloIdHitMap) {for (const auto& ptr : liste) seedList.insert(ptr);}
      if (diagLevel_ > 2) dump("Post filtering", caloIdHitMap,seedList);




      //produce split-offs clusters
      while (!seedList.empty())
      {
          const CaloHit* crystalSeed = *seedList.begin();
          ClusterFinder finder(cal,crystalSeed,deltaTime_, ExpandCut_, addSecondRing_);

          finder.formCluster(caloIdHitMap);
          splitClusterList.push_back(finder.clusterList());

          for (const auto& hit: finder.clusterList()) seedList.erase(hit);
      }

      //save the main and split clusters
      for (auto cluster : mainClusterList)  fillCluster(caloProtoClustersMain,cluster,CaloHitsHandle);
      for (auto cluster : splitClusterList) fillCluster(caloProtoClustersSplit,cluster,CaloHitsHandle);

      //sort these guys
      std::sort(caloProtoClustersMain.begin(),  caloProtoClustersMain.end(), [](const CaloProtoCluster& a, const CaloProtoCluster& b) {return a.time() < b.time();});
      std::sort(caloProtoClustersSplit.begin(), caloProtoClustersSplit.end(),[](const CaloProtoCluster& a, const CaloProtoCluster& b) {return a.time() < b.time();});

      int totNcrys(0); double totEnergy(0.0);
      for (const auto& clu : caloProtoClustersMain)  {totEnergy += clu.energyDep(); totNcrys +=clu.caloHitsPtrVector().size();}
      for (const auto& clu : caloProtoClustersSplit) {totEnergy += clu.energyDep(); totNcrys +=clu.caloHitsPtrVector().size();}
      if (diagLevel_ > 0) std::cout<<"[CaloProtoClusterMaker] Total energy / Ncrys ="<<totEnergy<<" "<<totNcrys<<std::endl;
  }





  //----------------------------------------------------------------------------------------------------------
  void CaloProtoClusterMaker::fillCluster(CaloProtoClusterCollection& caloProtoClustersColl,
                                          const CaloCrystalList& clusterPtrList,
                                          const art::Handle<CaloHitCollection>& CaloHitsHandle)
  {
      const CaloHitCollection& CaloHits(*CaloHitsHandle);
      const CaloHit* CaloHitBase = &CaloHits.front();

      std::vector<art::Ptr<CaloHit>> caloHitsPtrVector;
      double totalEnergy(0),totalEnergyErr(0);
      //double timeW(0),timeWtot(0);

      for (auto clusterPrt : clusterPtrList)
      {
          //double weight = 1.0/clusterPrt->timeErr()/clusterPrt->timeErr();
          //timeW    += weight*clusterPrt->time();
          //timeWtot += weight;

          totalEnergy    += clusterPrt->energyDep();
          totalEnergyErr += clusterPrt->energyDepErr()*clusterPrt->energyDepErr();

          size_t idx = (clusterPrt - CaloHitBase);
          caloHitsPtrVector.push_back(art::Ptr<CaloHit>(CaloHitsHandle,idx));
      }

      totalEnergyErr = sqrt(totalEnergyErr);
      double time    = (*clusterPtrList.begin())->time();
      double timeErr = (*clusterPtrList.begin())->timeErr();
      //double time    = timeW/timeWtot;
      //double timeErr = 1.0/sqrt(timeWtot);

      caloProtoClustersColl.emplace_back(CaloProtoCluster(time,timeErr,totalEnergy,totalEnergyErr,caloHitsPtrVector,false));

      if (diagLevel_ > 1)
      {
          std::cout<<"This cluster contains "<<clusterPtrList.size()<<" crystals, id= ";
          for (auto clusterPrt : clusterPtrList) std::cout<<clusterPrt->crystalID()<<" ";
          std::cout<<" with energy="<<totalEnergy<<" and time="<<time<<std::endl;;
      }
  }




  //----------------------------------------------------------------------------------------------------------
  void CaloProtoClusterMaker::filterByTime(CaloCrystalList& liste, const std::vector<double>& clusterTime)
  {
      auto it = liste.begin();
      while (it != liste.end())
      {
          const CaloHit* hit = *it;

          auto itTime = clusterTime.begin();
          while (itTime != clusterTime.end()){if ( (*itTime - hit->time()) < deltaTime_) break; ++itTime;}

          if (itTime == clusterTime.end() ) it = liste.erase(it);
          else ++it;
      }
  }



  //----------------------------------------------------------------------------------------------------------
  void CaloProtoClusterMaker::dump(const std::string& title, const std::vector<CaloCrystalList>& caloIdHitMap,
                                            std::set<const CaloHit*> seedList)
  {
      std::cout<<title<<std::endl;
      std::cout<<"Cache content"<<std::endl;
      for (unsigned i=0;i<caloIdHitMap.size();++i)
      {
         if (caloIdHitMap[i].empty()) continue;
         std::cout<<"Crystal idx "<<i<<std::endl;
         for (auto& ptr : caloIdHitMap[i]) std::cout<<ptr<<" "<<ptr->energyDep()<<"  ";
         std::cout<<std::endl;
      }
      std::cout<<"Seeds  "<<std::endl;
      for (auto& ptr :seedList) std::cout<<ptr<<" ";
      std::cout<<std::endl;
  }





}

DEFINE_ART_MODULE(mu2e::CaloProtoClusterMaker)
