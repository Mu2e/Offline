#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"

#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/Mu2eUtilities/inc/MVATools.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/TriggerInfo.hh"

#include <vector>
#include <string>



namespace mu2e {

  class FilterEcalNNTrigger : public art::EDFilter
  {
     public:
        struct Config
        {
            using Name    = fhicl::Name;
            using Comment = fhicl::Comment;
            fhicl::Atom<art::InputTag>     caloClusterCollection { Name("caloClusterCollection"), Comment("Calo cluster collection name") };
            fhicl::Table<MVATools::Config> caloBkgMVA            { Name("caloBkgMVA"),            Comment("MVA configuration") };
            fhicl::Atom<float>             minEtoTest            { Name("minEtoTest"),            Comment("Minimum energy to run the MVA") };
            fhicl::Atom<float>             minRtoTest            { Name("minRtoTest"),            Comment("Minimum radius to run the MVA") };
            fhicl::Atom<float>             minMVAScore           { Name("minMVAScore"),           Comment("MVA cut for signal") };
            fhicl::Atom<int>               diagLevel             { Name("diagLevel"),             Comment("Diag Level"),0 };
        };

        explicit FilterEcalNNTrigger(const art::EDFilter::Table<Config>& config) :
          EDFilter{config},
          caloClusterToken_{consumes<CaloClusterCollection>(config().caloClusterCollection())},
          caloBkgMVA_      (config().caloBkgMVA()),
          minEtoTest_      (config().minEtoTest()),
          minMVAScore_     (config().minMVAScore()),
          diagLevel_       (config().diagLevel())
        {
           produces<TriggerInfo>();
        }

        bool filter(art::Event& event) override;
        void beginJob() override;


     private:
        art::ProductToken<CaloClusterCollection> caloClusterToken_;
        MVATools caloBkgMVA_;
        float    minEtoTest_;
        float    minRtoTest_;
        float    minMVAScore_;
        int      diagLevel_;

        bool     filterClusters(const art::Handle<CaloClusterCollection>& caloClustersHandle, TriggerInfo& trigInfo);
        float    secondMoment(const Calorimeter& cal, const CaloHitPtrVector& hits) const;
  };


  void FilterEcalNNTrigger::beginJob()
  {
     caloBkgMVA_.initMVA();
  }



  bool FilterEcalNNTrigger::filter(art::Event& event)
  {
     art::Handle<CaloClusterCollection> caloClustersHandle = event.getHandle<CaloClusterCollection>(caloClusterToken_);

     auto trigInfo = std::make_unique<TriggerInfo>();
     bool retval   = filterClusters(caloClustersHandle, *trigInfo);
     event.put(std::move(trigInfo));

     return retval;
  }


  //----------------------------------------------------------------------------------------------------------
  bool FilterEcalNNTrigger::filterClusters(const art::Handle<CaloClusterCollection>& caloClustersHandle, TriggerInfo& trigInfo)
  {
     const Calorimeter& cal = *(GeomHandle<Calorimeter>());
     const CaloClusterCollection& caloClusters(*caloClustersHandle);

     bool select(false);
     for (auto clusterIt=caloClusters.begin(); clusterIt != caloClusters.end();++clusterIt){
        if (clusterIt->energyDep() < minEtoTest_)         continue;
        if (clusterIt->cog3Vector().perp() < minRtoTest_) continue;

        const auto& hits         = clusterIt->caloHitsPtrVector();
        const auto& neighborsId  = cal.crystal(hits[0]->crystalID()).neighbors();
        const auto& nneighborsId = cal.crystal(hits[0]->crystalID()).nextNeighbors();

        float e1(hits[0]->energyDep());
        float e2(hits[0]->energyDep());
        if (hits.size()>1) e2 += hits[1]->energyDep();

        float e9(hits[0]->energyDep()),e25(hits[0]->energyDep());
        for (auto hit : hits){
            if (std::find(neighborsId.begin(),  neighborsId.end(),  hit->crystalID()) != neighborsId.end())  {e9 += hit->energyDep();e25 += hit->energyDep();}
            if (std::find(nneighborsId.begin(), nneighborsId.end(), hit->crystalID()) != nneighborsId.end()) {e25 += hit->energyDep();}
        }

        float secondMom = secondMoment(cal,hits);

        std::vector<float> mvavars(8,0.0);
        //mvavars[0] = clusterIt->energyDep();
        mvavars[0] = clusterIt->cog3Vector().perp();
        mvavars[1] = clusterIt->time();
        mvavars[2] = clusterIt->size();
        mvavars[3] = e1/clusterIt->energyDep();
        mvavars[4] = e2/clusterIt->energyDep();
        mvavars[5] = e9/clusterIt->energyDep();
        mvavars[6] = e25/clusterIt->energyDep();
        mvavars[7] = secondMom;

        float mvaout = caloBkgMVA_.evalMVA(mvavars);
        if (mvaout < minMVAScore_) continue;

        select = true;
        size_t index = std::distance(caloClusters.begin(),clusterIt);
        trigInfo._caloClusters.push_back(art::Ptr<CaloCluster>(caloClustersHandle,index));
     }
     return select;
  }

  //-----------------------------------------------------------------------------------------------------------------------
  float FilterEcalNNTrigger::secondMoment(const Calorimeter& cal, const CaloHitPtrVector& hits) const
  {
     double sx(0),sy(0),sx2(0),sy2(0),sw(0);
     for (const auto& hit : hits){
        auto crId(hit->crystalID());
        auto energy(hit->energyDep());

        auto xCrystal = cal.crystal(crId).position().x();
        auto yCrystal = cal.crystal(crId).position().y();
        auto weight   = energy; //maybe try log(energy);

        sw  += weight;
        sx  += xCrystal*weight;
        sy  += yCrystal*weight;
        sx2 += xCrystal*xCrystal*weight;
        sy2 += yCrystal*yCrystal*weight;
     }
     return (sx2-sx*sx/sw + sy2-sy*sy/sw)/sw;
  }

}

DEFINE_ART_MODULE(mu2e::FilterEcalNNTrigger)
