#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"

#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/CaloCluster/inc/ClusterUtils.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/Mu2eUtilities/inc/MVATools.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/TriggerInfo.hh"
#include "Offline/RecoDataProducts/inc/MVAResult.hh"

#include <vector>
#include <string>


namespace mu2e {

  class CaloNNEval : public art::EDProducer
  {
     public:
        struct Config
        {
          using Name    = fhicl::Name;
          using Comment = fhicl::Comment;
          fhicl::Atom<art::InputTag>     caloClusterCollection { Name("caloClusterCollection"), Comment("Calo cluster collection name") };
          fhicl::Table<MVATools::Config> caloBkgMVA            { Name("caloBkgMVA"),            Comment("MVA configuration") };
          fhicl::Atom<float>             minEtoTest            { Name("minEtoTest"),            Comment("Minimum cluster energy to run the MVA") };
          fhicl::Atom<float>             minRtoTest            { Name("minRtoTest"),            Comment("Minimum cluster radius to run the MVA") };
          fhicl::Atom<float>             minTtoTest            { Name("minTtoTest"),            Comment("Minimum cluster time to run the MVA") };
          fhicl::Atom<float>             maxEtoTest            { Name("maxEtoTest"),            Comment("Maximum cluster energy to run the MVA") };
          fhicl::Atom<float>             maxRtoTest            { Name("maxRtoTest"),            Comment("Maximum cluster radius to run the MVA") };
          fhicl::Atom<float>             maxTtoTest            { Name("maxTtoTest"),            Comment("Maximum cluster time to run the MVA") };
          fhicl::Atom<int>               diagLevel             { Name("diagLevel"),             Comment("Diag Level"),0 };
        };

        explicit CaloNNEval(const  art::EDProducer::Table<Config>& config) :
          EDProducer{config},
          caloClusterToken_{consumes<CaloClusterCollection>(config().caloClusterCollection())},
          caloBkgMVA_      (config().caloBkgMVA()),
          minEtoTest_      (config().minEtoTest()),
          minRtoTest_      (config().minRtoTest()),
          minTtoTest_      (config().minTtoTest()),
          maxEtoTest_      (config().maxEtoTest()),
          maxRtoTest_      (config().maxRtoTest()),
          maxTtoTest_      (config().maxTtoTest()),
          diagLevel_       (config().diagLevel())
        {
          produces<MVAResultCollection>();
        }

        void produce(art::Event& event) override;
        void beginJob() override;


     private:
       art::ProductToken<CaloClusterCollection> caloClusterToken_;
       MVATools caloBkgMVA_;
       float    minEtoTest_;
       float    minRtoTest_;
       float    minTtoTest_;
       float    maxEtoTest_;
       float    maxRtoTest_;
       float    maxTtoTest_;
       int      diagLevel_;

       void     evalClusters(const art::Handle<CaloClusterCollection>& caloClustersHandle, MVAResultCollection& MVAresults);
       float    secondMoment(const Calorimeter& cal, const CaloHitPtrVector& hits) const;
  };


  void CaloNNEval::beginJob()
  {
    caloBkgMVA_.initMVA();
  }



  void CaloNNEval::produce(art::Event& event)
  {
    art::Handle<CaloClusterCollection> caloClustersHandle = event.getHandle<CaloClusterCollection>(caloClusterToken_);
    auto MVAresults = std::make_unique<MVAResultCollection>();
    evalClusters(caloClustersHandle, *MVAresults);
    event.put(std::move(MVAresults));
  }


  //----------------------------------------------------------------------------------------------------------
  void CaloNNEval::evalClusters(const art::Handle<CaloClusterCollection>& caloClustersHandle, MVAResultCollection& MVAresults)
  {
     if (!caloClustersHandle.isValid()) return;

     const Calorimeter& cal = *(GeomHandle<Calorimeter>());
     const CaloClusterCollection& caloClusters(*caloClustersHandle);

     for (auto clusterIt=caloClusters.begin(); clusterIt != caloClusters.end();++clusterIt){

        float cluE = clusterIt->energyDep();
        float cluR = clusterIt->cog3Vector().perp();
        float cluT = clusterIt->time();
        if (cluE < minEtoTest_ || cluE > maxEtoTest_ || cluR < minRtoTest_ || cluR > maxRtoTest_ ||
            cluT < minTtoTest_ || cluT > maxTtoTest_)
        {
          MVAresults.push_back(MVAResult(-1.0));
          if (diagLevel_>0) std::cout<<"CaloNNEval skip cluster with e="<<clusterIt->energyDep()
                                     <<" r="<< clusterIt->cog3Vector().perp()<<" t="<<clusterIt->time()<<"\n";
          continue;
        }

        ClusterUtils cluUtil(cal,*clusterIt);

        //This would be the version with energy and time included
        //std::vector<float> mvavars(9,0.0);
        //mvavars[0] = clusterIt->energyDep();
        //mvavars[0] = clusterIt->cog3Vector().perp();
        //mvavars[1] = clusterIt->time();
        //mvavars[2] = clusterIt->size();
        //mvavars[3] = cluUtil.e1()/clusterIt->energyDep();
        //mvavars[4] = cluUtil.e2()/clusterIt->energyDep();
        //mvavars[5] = cluUtil.e9()/clusterIt->energyDep();
        //mvavars[6] = cluUtil.e25()/clusterIt->energyDep();
        //mvavars[7] = cluUtil.secondMoment();

        std::vector<float> mvavars(7,0.0);
        mvavars[0] = clusterIt->cog3Vector().perp();
        mvavars[1] = clusterIt->size();
        mvavars[2] = cluUtil.e1()/clusterIt->energyDep();
        mvavars[3] = cluUtil.e2()/clusterIt->energyDep();
        mvavars[4] = cluUtil.e9()/clusterIt->energyDep();
        mvavars[5] = cluUtil.e25()/clusterIt->energyDep();
        mvavars[6] = cluUtil.secondMoment();

        float mvaout = caloBkgMVA_.evalMVA(mvavars);
        MVAresults.push_back(MVAResult(mvaout));
        if (diagLevel_>0) std::cout<<"CaloNNEval cluster with e="<<clusterIt->energyDep()
                                     <<" and r="<< clusterIt->cog3Vector().perp()
                                     <<"has mvaout="<<mvaout<<"\n";
     }
  }

}

DEFINE_ART_MODULE(mu2e::CaloNNEval)
