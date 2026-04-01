#include "art/Framework/Core/EDProducer.h"
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

        //This would be the version with energy and time included
        //std::vector<float> mvavars(9,0.0);
        //mvavars[0] = clusterIt->energyDep();
        //mvavars[0] = clusterIt->cog3Vector().perp();
        //mvavars[1] = clusterIt->time();
        //mvavars[2] = clusterIt->size();
        //mvavars[3] = e1/clusterIt->energyDep();
        //mvavars[4] = e2/clusterIt->energyDep();
        //mvavars[5] = e9/clusterIt->energyDep();
        //mvavars[6] = e25/clusterIt->energyDep();
        //mvavars[7] = secondMom;

        std::vector<float> mvavars(7,0.0);
        mvavars[0] = clusterIt->cog3Vector().perp();
        mvavars[1] = clusterIt->size();
        mvavars[2] = e1/clusterIt->energyDep();
        mvavars[3] = e2/clusterIt->energyDep();
        mvavars[4] = e9/clusterIt->energyDep();
        mvavars[5] = e25/clusterIt->energyDep();
        mvavars[6] = secondMom;

        float mvaout = caloBkgMVA_.evalMVA(mvavars);
        MVAresults.push_back(MVAResult(mvaout));
        if (diagLevel_>0) std::cout<<"CaloNNEval cluster with e="<<clusterIt->energyDep()
                                     <<" and r="<< clusterIt->cog3Vector().perp()
                                     <<"has mvaout="<<mvaout<<"\n";
     }
  }

  //-----------------------------------------------------------------------------------------------------------------------
  float CaloNNEval::secondMoment(const Calorimeter& cal, const CaloHitPtrVector& hits) const
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

DEFINE_ART_MODULE(mu2e::CaloNNEval)
