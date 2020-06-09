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
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"

#include "CaloCluster/inc/ClusterFinder.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "RecoDataProducts/inc/CaloCrystalHit.hh"
#include "RecoDataProducts/inc/CaloProtoCluster.hh"

#include <iostream>
#include <string>
#include <list>
#include <set>
#include <vector>


namespace mu2e {

  class CaloProtoClusterFromCrystalHit : public art::EDProducer
  {
     public:
        typedef std::vector<const CaloCrystalHit*>  CaloCrystalVec;
        typedef std::list<const CaloCrystalHit*>    CaloCrystalList;
        
        struct Config
        {
            using Name    = fhicl::Name;
            using Comment = fhicl::Comment;
            fhicl::Atom<std::string>  caloCrystalModuleLabel{ Name("caloCrystalModuleLabel"), Comment("Calo Crystal module label")};
            fhicl::Atom<std::string>  mainTag               { Name("mainTag"),                Comment("Name of main cluster collection") }; 
            fhicl::Atom<std::string>  splitTag              { Name("splitTag"),               Comment("NAme of split-off collection") }; 
            fhicl::Atom<double>       EminSeed              { Name("EminSeed"),               Comment("Minimum energy for a hit to be a cluster seed") }; 
            fhicl::Atom<double>       EnoiseCut             { Name("EnoiseCut"),              Comment("Minimum energy for a hit to be in a cluster") }; 
            fhicl::Atom<double>       ExpandCut             { Name("ExpandCut"),              Comment("Minimum energy for a hit to expand cluster") }; 
            fhicl::Atom<double>       timeCut               { Name("timeCut"),                Comment("Minimum hit time to form cluster") }; 
            fhicl::Atom<double>       deltaTime             { Name("deltaTime"),              Comment("Maximum time difference between seed and hit in cluster") }; 
            fhicl::Atom<int>          diagLevel             { Name("diagLevel"),              Comment("Diag level"),0 }; 
        };

        explicit CaloProtoClusterFromCrystalHit(const art::EDProducer::Table<Config>& config) :
          EDProducer{config},
          caloCrystalToken_{consumes<CaloCrystalHitCollection>(config().caloCrystalModuleLabel())},
          mainTag_         (config().mainTag()),
          splitTag_        (config().splitTag()),
          EminSeed_        (config().EminSeed()),
          EnoiseCut_       (config().EnoiseCut()),
          ExpandCut_       (config().ExpandCut()),
          timeCut_         (config().timeCut()),
          deltaTime_       (config().deltaTime()),
          diagLevel_       (config().diagLevel())
        {
           produces<CaloProtoClusterCollection>(mainTag_);
           produces<CaloProtoClusterCollection>(splitTag_);
        }

        void produce(art::Event& e) override;

     private:
        art::ProductToken<CaloCrystalHitCollection> caloCrystalToken_;
        std::string       mainTag_;
        std::string       splitTag_;
        double            EminSeed_;
        double            EnoiseCut_;
        double            ExpandCut_;
        double            timeCut_;
        double            deltaTime_;
        int               diagLevel_;

        void makeProtoClusters (CaloProtoClusterCollection&,CaloProtoClusterCollection&, const art::Handle<CaloCrystalHitCollection>&);
        void filterByTime      (CaloCrystalList&, const std::vector<double>&);
        void fillCluster       (CaloProtoClusterCollection&, const CaloCrystalList&,const art::Handle<CaloCrystalHitCollection>&);
        void dump              (const std::string&, const std::vector<CaloCrystalList>&, std::set<const CaloCrystalHit*>);
  };


  void CaloProtoClusterFromCrystalHit::produce(art::Event& event)
  {
      // Check that calorimeter geometry description exists
      art::ServiceHandle<GeometryService> geom;
      if( !(geom->hasElement<Calorimeter>()) ) return;

      // Get handles to calorimeter crystal hits
      art::Handle<CaloCrystalHitCollection> CaloCrystalHitsHandle;
      bool const success = event.getByToken(caloCrystalToken_, CaloCrystalHitsHandle);
      if (!success) return;

      // Create a new CaloCluster collection and fill it
      auto caloProtoClustersMain = std::make_unique<CaloProtoClusterCollection>();
      auto caloProtoClustersSplit = std::make_unique<CaloProtoClusterCollection>();
      makeProtoClusters(*caloProtoClustersMain,*caloProtoClustersSplit,CaloCrystalHitsHandle);

      event.put(std::move(caloProtoClustersMain),  mainTag_);
      event.put(std::move(caloProtoClustersSplit), splitTag_);
  }


  //----------------------------------------------------------------------------------------------------------
  void CaloProtoClusterFromCrystalHit::makeProtoClusters(CaloProtoClusterCollection& caloProtoClustersMain,
                                                         CaloProtoClusterCollection& caloProtoClustersSplit,
                                                         const art::Handle<CaloCrystalHitCollection> & CaloCrystalHitsHandle)
  {
      const Calorimeter& cal = *(GeomHandle<Calorimeter>());
      const CaloCrystalHitCollection& CaloCrystalHits(*CaloCrystalHitsHandle);
      if (CaloCrystalHits.empty()) return;


      //declare and fill the hash map crystal_id -> list of CaloHits
      std::vector<CaloCrystalList>    mainClusterList, splitClusterList, caloIdHitMap(cal.nCrystal());
      std::set<const CaloCrystalHit*> seedList;
      std::vector<double>             clusterTime;

      //fill data structures
      for (const auto& hit : CaloCrystalHits)
      {
          if (hit.energyDep() < EnoiseCut_ || hit.time() < timeCut_) continue;
          caloIdHitMap[hit.id()].push_back(&hit);
          if (hit.energyDep() > EminSeed_ ) seedList.insert(&hit);
      }
      
      if (diagLevel_ > 2) dump("Init", caloIdHitMap,seedList);
       


      //produce main clusters
      while( !seedList.empty() )
      {
          const CaloCrystalHit* crystalSeed = *seedList.begin();
          ClusterFinder finder(cal, crystalSeed, deltaTime_, ExpandCut_);
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
          const CaloCrystalHit* crystalSeed = *seedList.begin();
          ClusterFinder finder(cal,crystalSeed,deltaTime_, ExpandCut_);

          finder.formCluster(caloIdHitMap);
          splitClusterList.push_back(finder.clusterList());

          for (const auto& hit: finder.clusterList()) seedList.erase(hit);
      }

      //save the main and split clusters
      for (auto cluster : mainClusterList)  fillCluster(caloProtoClustersMain,cluster,CaloCrystalHitsHandle);
      for (auto cluster : splitClusterList) fillCluster(caloProtoClustersSplit,cluster,CaloCrystalHitsHandle);

      //sort these guys
      std::sort(caloProtoClustersMain.begin(),  caloProtoClustersMain.end(), [](const CaloProtoCluster& a, const CaloProtoCluster& b) {return a.time() < b.time();});
      std::sort(caloProtoClustersSplit.begin(), caloProtoClustersSplit.end(),[](const CaloProtoCluster& a, const CaloProtoCluster& b) {return a.time() < b.time();});

      int totNcrys(0); double totEnergy(0.0);
      for (const auto& clu : caloProtoClustersMain)  {totEnergy += clu.energyDep(); totNcrys +=clu.caloCrystalHitsPtrVector().size();}
      for (const auto& clu : caloProtoClustersSplit) {totEnergy += clu.energyDep(); totNcrys +=clu.caloCrystalHitsPtrVector().size();}
      if (diagLevel_ > 0) std::cout<<"[CaloProtoClusterFromCrystalHit] Total energy / Ncrys ="<<totEnergy<<" "<<totNcrys<<std::endl;
  }





  //----------------------------------------------------------------------------------------------------------
  void CaloProtoClusterFromCrystalHit::fillCluster(CaloProtoClusterCollection& caloProtoClustersColl, 
                                                   const CaloCrystalList& clusterPtrList,
                                                   const art::Handle<CaloCrystalHitCollection>& CaloCrystalHitsHandle)
  {
      const CaloCrystalHitCollection& CaloCrystalHits(*CaloCrystalHitsHandle);
      const CaloCrystalHit* caloCrystalHitBase = &CaloCrystalHits.front();

      std::vector<art::Ptr<CaloCrystalHit>> caloCrystalHitsPtrVector;
      double totalEnergy(0),totalEnergyErr(0);
      //double timeW(0),timeWtot(0);

      for (auto clusterPrt : clusterPtrList)
      {
          //double weight = 1.0/clusterPrt->timeErr()/clusterPrt->timeErr();
          //timeW    += weight*clusterPrt->time();
          //timeWtot += weight;

          totalEnergy    += clusterPrt->energyDep();
          totalEnergyErr += clusterPrt->energyDepErr()*clusterPrt->energyDepErr();

          size_t idx = (clusterPrt - caloCrystalHitBase);
          caloCrystalHitsPtrVector.push_back(art::Ptr<CaloCrystalHit>(CaloCrystalHitsHandle,idx));
      }

      totalEnergyErr = sqrt(totalEnergyErr);
      double time    = (*clusterPtrList.begin())->time();
      double timeErr = (*clusterPtrList.begin())->timeErr();
      //double time    = timeW/timeWtot;
      //double timeErr = 1.0/sqrt(timeWtot);

      caloProtoClustersColl.emplace_back(CaloProtoCluster(time,timeErr,totalEnergy,totalEnergyErr,caloCrystalHitsPtrVector,false));

      if (diagLevel_ > 1)
      {
          std::cout<<"This cluster contains "<<clusterPtrList.size()<<" crystals, id= ";
          for (auto clusterPrt : clusterPtrList) std::cout<<clusterPrt->id()<<" ";
          std::cout<<" with energy="<<totalEnergy<<" and time="<<time<<std::endl;;
      }
  }




  //----------------------------------------------------------------------------------------------------------
  void CaloProtoClusterFromCrystalHit::filterByTime(CaloCrystalList& liste, const std::vector<double>& clusterTime)
  {
      auto it = liste.begin();
      while (it != liste.end())
      {
          const CaloCrystalHit* hit = *it;
          
          auto itTime = clusterTime.begin();
          while (itTime != clusterTime.end()){if ( (*itTime - hit->time()) < deltaTime_) break; ++itTime;}
                    
          if (itTime == clusterTime.end() ) it = liste.erase(it);
          else ++it;
      }
  }



  //----------------------------------------------------------------------------------------------------------
  void CaloProtoClusterFromCrystalHit::dump(const std::string& title, const std::vector<CaloCrystalList>& caloIdHitMap, 
                                            std::set<const CaloCrystalHit*> seedList)
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

DEFINE_ART_MODULE(mu2e::CaloProtoClusterFromCrystalHit);
