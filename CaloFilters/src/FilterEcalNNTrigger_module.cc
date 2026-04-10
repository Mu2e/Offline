#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"

#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/MVAResult.hh"
#include "Offline/RecoDataProducts/inc/TriggerInfo.hh"


namespace mu2e {

  class FilterEcalNNTrigger : public art::EDFilter
  {
     public:
        struct Config
        {
            using Name    = fhicl::Name;
            using Comment = fhicl::Comment;
            fhicl::Atom<art::InputTag>     caloClusterCollection { Name("caloClusterCollection"), Comment("Calo cluster collection name") };
            fhicl::Atom<art::InputTag>     caloMVACollection     { Name("caloMVACollection"),     Comment("Calo MVA cluster collection name") };
            fhicl::Atom<float>             minEtoTest            { Name("minEtoTest"),            Comment("Minimum cluster energy to run the MVA") };
            fhicl::Atom<float>             minRtoTest            { Name("minRtoTest"),            Comment("Minimum cluster radius to run the MVA") };
            fhicl::Atom<float>             minTtoTest            { Name("minTtoTest"),            Comment("Minimum cluster time to run the MVA") };
            fhicl::Atom<float>             maxEtoTest            { Name("maxEtoTest"),            Comment("Maximum cluster energy to run the MVA") };
            fhicl::Atom<float>             maxRtoTest            { Name("maxRtoTest"),            Comment("Maximum cluster radius to run the MVA") };
            fhicl::Atom<float>             maxTtoTest            { Name("maxTtoTest"),            Comment("Maximum cluster time to run the MVA") };
            fhicl::Atom<float>             minMVAScore           { Name("minMVAScore"),           Comment("MVA cut for signal") };
            fhicl::Atom<int>               diagLevel             { Name("diagLevel"),             Comment("Diag Level"),0 };
        };

        explicit FilterEcalNNTrigger(const art::EDFilter::Table<Config>& config) :
          EDFilter{config},
          caloClusterToken_{consumes<CaloClusterCollection>(config().caloClusterCollection())},
          caloMVAToken_    {consumes<MVAResultCollection>(config().caloMVACollection())},
          minEtoTest_      (config().minEtoTest()),
          minRtoTest_      (config().minRtoTest()),
          minTtoTest_      (config().minTtoTest()),
          maxEtoTest_      (config().maxEtoTest()),
          maxRtoTest_      (config().maxRtoTest()),
          maxTtoTest_      (config().maxTtoTest()),
          minMVAScore_     (config().minMVAScore()),
          diagLevel_       (config().diagLevel())
        {
           produces<TriggerInfo>();
        }

        bool filter(art::Event& event) override;


     private:
        art::ProductToken<CaloClusterCollection> caloClusterToken_;
        art::ProductToken<MVAResultCollection>   caloMVAToken_;
        float    minEtoTest_;
        float    minRtoTest_;
        float    minTtoTest_;
        float    maxEtoTest_;
        float    maxRtoTest_;
        float    maxTtoTest_;
        float    minMVAScore_;
        int      diagLevel_;

        bool     filterClusters(const art::Handle<CaloClusterCollection>& caloClustersHandle,
                                const art::Handle<MVAResultCollection>& caloMVAHandle,
                                TriggerInfo& trigInfo);
  };


  bool FilterEcalNNTrigger::filter(art::Event& event)
  {
     art::Handle<CaloClusterCollection> caloClustersHandle = event.getHandle<CaloClusterCollection>(caloClusterToken_);
     art::Handle<MVAResultCollection>   caloMVAHandle      = event.getHandle<MVAResultCollection>(caloMVAToken_);

     auto trigInfo = std::make_unique<TriggerInfo>();
     bool retval   = filterClusters(caloClustersHandle, caloMVAHandle, *trigInfo);
     event.put(std::move(trigInfo));

     return retval;
  }


  //----------------------------------------------------------------------------------------------------------
  bool FilterEcalNNTrigger::filterClusters(const art::Handle<CaloClusterCollection>& caloClustersHandle,
                                           const art::Handle<MVAResultCollection>& caloMVAHandle,
                                           TriggerInfo& trigInfo)
  {
     if (!caloClustersHandle.isValid() || !caloMVAHandle.isValid()) return false;

     const auto& caloClusters(*caloClustersHandle);
     const auto& caloMVAs(*caloMVAHandle);

     if (caloClusters.size() != caloMVAs.size()){
       throw cet::exception("FILTER")<< "FilterEcalNNTrigger: Clusters and MVA collection sizes incpmpatible\n";
     }

     bool select(false);
     for (size_t index=0;index<caloClusters.size();++index){
        const auto& cluster = caloClusters[index];
        const auto& mvaout  = caloMVAs[index];

        float cluE = cluster.energyDep();
        float cluR = cluster.cog3Vector().perp();
        float cluT = cluster.time();
        if (cluE < minEtoTest_ || cluE > maxEtoTest_) continue;
        if (cluR < minRtoTest_ || cluR > maxRtoTest_) continue;
        if (cluT < minTtoTest_ || cluT > maxTtoTest_) continue;
        if (mvaout._value < minMVAScore_)             continue;

        select = true;
        trigInfo._caloClusters.push_back(art::Ptr<CaloCluster>(caloClustersHandle,index));
     }
     return select;
  }

}

DEFINE_ART_MODULE(mu2e::FilterEcalNNTrigger)
