//
// Module to produce the calorimeter clusters from proto-clusters.
//
// The strategy is to attach the split-off to the main cluster first, take the closest cluster if the split-off time is compatible
// with several clusters. Then associate energetic clusters between them, including the reattached split-off of each cluster in the
// comparison.
//
// Note: the cluster center-of-gravity is calculated in the calorimeter section front face frame
//
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"

#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/CaloCluster/inc/ClusterAssociator.hh"
#include "Offline/CaloCluster/inc/ClusterUtils.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/CaloProtoCluster.hh"

#include <iostream>
#include <string>
#include <vector>
#include <array>



namespace mu2e {

  class CaloClusterMaker : public art::EDProducer
  {
     public:
        struct Config
        {
            using Name    = fhicl::Name;
            using Comment = fhicl::Comment;
            fhicl::Atom<art::InputTag> caloClusterMainCollection  { Name("caloClusterMainCollection"),  Comment("CaloProtoClustere main collection name ") };
            fhicl::Atom<art::InputTag> caloClusterSplitCollection { Name("caloClusterSplitCollection"), Comment("CaloProtoClustere main collection name ") };
            fhicl::Atom<double>        deltaTime                  { Name("deltaTime"),                  Comment("Maximum time difference to associate clusters") };
            fhicl::Atom<double>        maxDistSplit               { Name("maxDistSplit"),               Comment("Maximum distance between split off and main cluster") };
            fhicl::Atom<double>        maxDistMain                { Name("maxDistMain"),                Comment("Maximum distance between main clusters") };
            fhicl::Atom<int>           strategy                   { Name("strategy"),                   Comment("Main cluster associator strategy") };
            fhicl::Atom<int>           diagLevel                  { Name("diagLevel"),                  Comment("Diag level"),0 };
        };

        explicit CaloClusterMaker(const art::EDProducer::Table<Config>& config) :
          EDProducer{config},
          mainToken_    {consumes<CaloProtoClusterCollection>(config().caloClusterMainCollection())},
          splitToken_   {consumes<CaloProtoClusterCollection>(config().caloClusterSplitCollection())},
          deltaTime_    (config().deltaTime()),
          maxDistSplit_ (config().maxDistSplit()),
          maxDistMain_  (config().maxDistMain()),
          strategy_     (config().strategy()),
          diagLevel_    (config().diagLevel())
        {
            produces<CaloClusterCollection>();
        }

        void produce(art::Event& e) override;


     private:
        art::ProductToken<CaloProtoClusterCollection>  mainToken_;
        art::ProductToken<CaloProtoClusterCollection>  splitToken_;
        double                                         deltaTime_;
        double                                         maxDistSplit_;
        double                                         maxDistMain_;
        int                                            strategy_;
        int                                            diagLevel_;

        std::vector<double> clusterTimeEnergy(const std::vector<art::Ptr<CaloHit>>& hits);
        void                makeCaloClusters(CaloClusterCollection&, const CaloProtoClusterCollection&, const CaloProtoClusterCollection&);
  };


  //---------------------------------------------------------------------------------------------------------------
  void CaloClusterMaker::produce(art::Event& event)
  {
      const auto& caloClustersMain  = *event.getValidHandle(mainToken_);
      const auto& caloClustersSplit = *event.getValidHandle(splitToken_);

      auto caloClusters = std::make_unique<CaloClusterCollection>();
      makeCaloClusters(*caloClusters, caloClustersMain, caloClustersSplit);

      event.put(std::move(caloClusters));
  }


  //---------------------------------------------------------------------------------------------------------------
  void CaloClusterMaker::makeCaloClusters(CaloClusterCollection& caloClusters, const CaloProtoClusterCollection& caloClustersMain,
                                          const CaloProtoClusterCollection& caloClustersSplit)
  {
      const Calorimeter& cal = *(GeomHandle<Calorimeter>());


      ClusterAssociator associator(cal);
      associator.associateSplitOff(caloClustersMain, caloClustersSplit, deltaTime_,maxDistSplit_);


      //-- First, associate split-off to main cluster into intermediate buffer
      CaloProtoClusterCollection caloProtoClustersTemp;

      for (unsigned imain(0); imain<caloClustersMain.size(); ++imain)
      {
          bool isSplit(false);
          std::vector<art::Ptr<CaloHit>> caloHitsPtrVector = caloClustersMain.at(imain).caloHitsPtrVector();

          //search split-offs and add their hits into the main cluster
          for (unsigned isplit=0;isplit<caloClustersSplit.size();++isplit)
          {
              if (associator.associatedSplitId(isplit) != imain) continue;
              isSplit = true;

              caloHitsPtrVector.insert(caloHitsPtrVector.end(),
                                              caloClustersSplit.at(isplit).caloHitsPtrVector().begin(),
                                              caloClustersSplit.at(isplit).caloHitsPtrVector().end());

              if (diagLevel_ > 1) std::cout<<"Associated main cluster "<<imain <<" with split cluster "<<isplit<<std::endl;
          }

          auto timeEnergy = clusterTimeEnergy(caloHitsPtrVector);
          caloProtoClustersTemp.emplace_back(CaloProtoCluster(timeEnergy[0],timeEnergy[1],timeEnergy[2],timeEnergy[3],caloHitsPtrVector,isSplit));
      }



      // combine main clusters together (split-off included in main clusters at this point)
      associator.associateMain(caloProtoClustersTemp, deltaTime_, maxDistMain_, strategy_);


      //finally, form final clusters
      std::vector<int> flagProto(caloClustersMain.size(),0);
      for (unsigned iproto=0;iproto<caloProtoClustersTemp.size();++iproto)
      {
          if (flagProto[iproto]) continue;

          auto        caloHitsPtrVector = caloProtoClustersTemp.at(iproto).caloHitsPtrVector();
          bool        isSplit           = caloProtoClustersTemp.at(iproto).isSplit();
          const auto& seed              = **caloHitsPtrVector.begin();
          int         diskId            = cal.crystal(seed.crystalID()).diskID();

          for (int iassoc : associator.associatedMainId(iproto))
          {
              flagProto[iassoc] = 1;
              isSplit           = true;

              caloHitsPtrVector.insert(caloHitsPtrVector.end(),
                                       caloProtoClustersTemp.at(iassoc).caloHitsPtrVector().begin(),
                                       caloProtoClustersTemp.at(iassoc).caloHitsPtrVector().end());

              if (diagLevel_ > 1) std::cout<<"Associated to main cluster id="<<iproto<<"   main split="<<iassoc<<std::endl;
          }

          //sort the crystal by energy
          auto cmpEnergy = [](const art::Ptr<CaloHit>& lhs,const art::Ptr<CaloHit>& rhs) {return lhs->energyDep() > rhs->energyDep();};
          std::sort(caloHitsPtrVector.begin(),caloHitsPtrVector.end(),cmpEnergy);


          ClusterUtils cluUtil(cal,caloHitsPtrVector);
          const auto cog = cluUtil.cog3Vector();
          auto timeEnergy = clusterTimeEnergy(caloHitsPtrVector);

          CaloCluster caloCluster(diskId,timeEnergy[0],timeEnergy[1],timeEnergy[2],timeEnergy[3],cog, caloHitsPtrVector,caloHitsPtrVector.size(),isSplit);
          caloClusters.push_back(caloCluster);

          if (diagLevel_ > 2)
          {
              std::cout<<"Making a new cluster with id= ";
              for (auto il = caloHitsPtrVector.begin(); il !=caloHitsPtrVector.end(); ++il) std::cout<<(*il)->crystalID()<<" ";
              std::cout<<std::endl;
          }
      }

      //finally, sort clusters by energy
      auto cmpCluster = [](CaloCluster& lhs, CaloCluster& rhs) {return lhs.energyDep() > rhs.energyDep();};
      std::sort(caloClusters.begin(),caloClusters.end(),cmpCluster);
  }


  //-------------------------------------------------------------------------------------------------
  std::vector<double> CaloClusterMaker::clusterTimeEnergy(const std::vector<art::Ptr<CaloHit>>& hits)
  {
      double totalEnergy(0),totalEnergyErr(0);

      for (auto hit : hits)
      {
          totalEnergy += hit->energyDep();
          totalEnergyErr += hit->energyDepErr()*hit->energyDepErr();
      }

      totalEnergyErr = sqrt(totalEnergyErr);

      const auto& seed = *hits.begin();
      double time      = seed->time();
      double timeErr   = seed->timeErr();

      return std::vector<double>{time,timeErr,totalEnergy,totalEnergyErr};
  }


}

DEFINE_ART_MODULE(mu2e::CaloClusterMaker)





/*
  To calculate clsuter time with weighted mean

  for (auto hit : hits)
  {
     double weight = 1.0/hit->timeErr()/hit->timeErr();
     timeW    += weight*hit->time();
     timeWtot += weight;
  }
  double time = timeW/timeWtot;
  double timeErr = 1.0/sqrt(timeWtot);
*/
