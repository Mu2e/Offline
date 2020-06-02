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
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"

#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CaloCluster/inc/ClusterAssociator.hh"
#include "CaloCluster/inc/ClusterMoments.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "RecoDataProducts/inc/CaloProtoClusterCollection.hh"

#include <iostream>
#include <string>
#include <vector>
#include <array>



namespace mu2e {

  class CaloClusterFromProtoCluster : public art::EDProducer 
  {
     public:
        struct Config
        {
            using Name    = fhicl::Name;
            using Comment = fhicl::Comment;
            fhicl::Atom<std::string>  caloClusterModuleLabel { Name("caloClusterModuleLabel"), Comment("Calo proto cluster collection name ") }; 
            fhicl::Atom<std::string>  mainTag                { Name("mainTag"),                Comment("Main cluster colletion tag name") }; 
            fhicl::Atom<std::string>  splitTag               { Name("splitTag"),               Comment("SPlit-off cluster colletion tag name") }; 
            fhicl::Atom<double>       deltaTime              { Name("deltaTime"),              Comment("Maximum time difference to associate clusters") }; 
            fhicl::Atom<double>       maxDistSplit           { Name("maxDistSplit"),           Comment("Maximum distance between split off and main cluster") }; 
            fhicl::Atom<double>       maxDistMain            { Name("maxDistMain"),            Comment("Maximum distance between main clusters") }; 
            fhicl::Atom<std::string>  cogTypeName            { Name("cogTypeName"),            Comment("Center-of-gravity calculator method") }; 
            fhicl::Atom<int>          strategy               { Name("strategy"),               Comment("Main cluster associator strategy") }; 
            fhicl::Atom<int>          diagLevel              { Name("diagLevel"),              Comment("Diag level "),0 }; 
        };
        
        explicit CaloClusterFromProtoCluster(const art::EDProducer::Table<Config>& config) :
          EDProducer{config},
          mainTag_      {config().caloClusterModuleLabel(), config().mainTag()},
          splitTag_     {config().caloClusterModuleLabel(), config().splitTag()},
          mainToken_    {consumes<CaloProtoClusterCollection>(mainTag_)},
          splitToken_   {consumes<CaloProtoClusterCollection>(splitTag_)},
          deltaTime_    (config().deltaTime()),
          maxDistSplit_ (config().maxDistSplit()),
          maxDistMain_  (config().maxDistMain()),
          cogTypeName_  (config().cogTypeName()),
          strategy_     (config().strategy()),
          diagLevel_    (config().diagLevel())
        {
            produces<CaloClusterCollection>();
        }

        void beginJob() override;
        void produce(art::Event& e) override;


     private:
        art::InputTag                                  mainTag_;
        art::InputTag                                  splitTag_;
        art::ProductToken<CaloProtoClusterCollection>  mainToken_;
        art::ProductToken<CaloProtoClusterCollection>  splitToken_;
        double                                         deltaTime_;
        double                                         maxDistSplit_;
        double                                         maxDistMain_;
        std::string                                    cogTypeName_;
        ClusterMoments::cogtype                        cogType_;
        int                                            strategy_;
        int                                            diagLevel_;

        void                 makeCaloClusters  (CaloClusterCollection&, const CaloProtoClusterCollection&, const CaloProtoClusterCollection&);
        std::array<double,3> calcEnergyLayer   (const Calorimeter&,const std::vector<art::Ptr<CaloCrystalHit>>&);
        std::array<double,4> clusterTimeEnergy (const std::vector<art::Ptr<CaloCrystalHit>>& hits);
  };



  //---------------------------------------------------------------------------------------------------------------
  void CaloClusterFromProtoCluster::beginJob()
  {
      if (cogTypeName_.compare("Linear"))    cogType_ = ClusterMoments::Linear;
      if (cogTypeName_.compare("Logarithm")) cogType_ = ClusterMoments::Logarithm;
  }

  
  //---------------------------------------------------------------------------------------------------------------
  void CaloClusterFromProtoCluster::produce(art::Event& event)
  {
      // Check that calorimeter geometry description exists
      art::ServiceHandle<GeometryService> geom;
      if( !(geom->hasElement<Calorimeter>()) ) return;

      const auto& caloClustersMain = *event.getValidHandle(mainToken_);
      const auto& caloClustersSplit = *event.getValidHandle(splitToken_);

      //Create a new CaloCluster collection and fill it
      auto caloClusters = std::make_unique<CaloClusterCollection>();
      makeCaloClusters(*caloClusters, caloClustersMain, caloClustersSplit);

      event.put(std::move(caloClusters));
  }


  //---------------------------------------------------------------------------------------------------------------
  void CaloClusterFromProtoCluster::makeCaloClusters(CaloClusterCollection& caloClusters,
                                                     const CaloProtoClusterCollection& caloClustersMain,
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
          std::vector<art::Ptr<CaloCrystalHit>> caloCrystalHitsPtrVector = caloClustersMain.at(imain).caloCrystalHitsPtrVector();

          //search split-offs and add their hits into the main cluster
          for (unsigned isplit=0;isplit<caloClustersSplit.size();++isplit)
          {
              if (associator.associatedSplitId(isplit) != imain) continue;
              isSplit = true;

              caloCrystalHitsPtrVector.insert(caloCrystalHitsPtrVector.end(),
                                              caloClustersSplit.at(isplit).caloCrystalHitsPtrVector().begin(),
                                              caloClustersSplit.at(isplit).caloCrystalHitsPtrVector().end());

              if (diagLevel_ > 1) std::cout<<"Associated main cluster "<<imain <<" with split cluster "<<isplit<<std::endl;
          }

          auto timeEnergy = clusterTimeEnergy(caloCrystalHitsPtrVector);
          caloProtoClustersTemp.emplace_back(CaloProtoCluster(timeEnergy[0],timeEnergy[1],timeEnergy[2],timeEnergy[3],caloCrystalHitsPtrVector,isSplit));
      }



      // combine main clusters together (split-off included in main clusters at this point)
      associator.associateMain(caloProtoClustersTemp, deltaTime_, maxDistMain_, strategy_);




      //finally, form final clusters
      std::vector<int> flagProto(caloClustersMain.size(),0);
      for (unsigned iproto=0;iproto<caloProtoClustersTemp.size();++iproto)
      {
          if (flagProto[iproto]) continue;

          auto        caloCrystalHitsPtrVector = caloProtoClustersTemp.at(iproto).caloCrystalHitsPtrVector();
          bool        isSplit                  = caloProtoClustersTemp.at(iproto).isSplit();
          const auto& seed                     = **caloCrystalHitsPtrVector.begin();
          int         diskId                   = cal.crystal(seed.id()).diskId();

          for (int iassoc : associator.associatedMainId(iproto))
          {
              flagProto[iassoc] = 1;
              isSplit           = true;

              caloCrystalHitsPtrVector.insert(caloCrystalHitsPtrVector.end(),
                                              caloProtoClustersTemp.at(iassoc).caloCrystalHitsPtrVector().begin(),
                                              caloProtoClustersTemp.at(iassoc).caloCrystalHitsPtrVector().end());

              if (diagLevel_ > 1) std::cout<<"Associated to main cluster id="<<iproto<<"   main split="<<iassoc<<std::endl;
          }

          //sort the crystal by energy
          auto cmpEnergy = [](const art::Ptr<CaloCrystalHit>& lhs,const art::Ptr<CaloCrystalHit>& rhs) {return lhs->energyDep() > rhs->energyDep();};
          std::sort(caloCrystalHitsPtrVector.begin(),caloCrystalHitsPtrVector.end(),cmpEnergy);


          auto timeEnergy = clusterTimeEnergy(caloCrystalHitsPtrVector);
          CaloCluster caloCluster(diskId,timeEnergy[0],timeEnergy[1],timeEnergy[2],timeEnergy[3],caloCrystalHitsPtrVector,caloCrystalHitsPtrVector.size(),isSplit);

          //calculate a lot of fancy useful things
          auto EnerLayer = calcEnergyLayer(cal,caloCrystalHitsPtrVector);
          ClusterMoments cogCalculator(cal,caloCluster,diskId);
          cogCalculator.calculate(cogType_);
          caloCluster.cog3Vector(cogCalculator.cog());
          caloCluster.secondMoment(cogCalculator.secondMoment());
          caloCluster.angle(cogCalculator.angle());
          caloCluster.energyRing(EnerLayer[0],EnerLayer[1],EnerLayer[2]);


          caloClusters.push_back(caloCluster);

          if (diagLevel_ > 2)
          {
              std::cout<<"Making a new cluster with id= ";
              for (auto il = caloCrystalHitsPtrVector.begin(); il !=caloCrystalHitsPtrVector.end(); ++il) std::cout<<(*il)->id()<<" ";
              std::cout<<std::endl;
          }
      }

      //finally, sort clusters by energy
      auto cmpCluster = [](CaloCluster& lhs, CaloCluster& rhs) {return lhs.energyDep() > rhs.energyDep();};
      std::sort(caloClusters.begin(),caloClusters.end(),cmpCluster);
  }




  //----------------------------------------------------------------------------------------------------
  std::array<double,3> CaloClusterFromProtoCluster::calcEnergyLayer(const Calorimeter& cal, const std::vector<art::Ptr<CaloCrystalHit>>& caloCrystalHitsPtrVector)
  {
      int seedId              = caloCrystalHitsPtrVector[0]->id();
      double seedEnergy       = caloCrystalHitsPtrVector[0]->energyDep();
      const auto neighborsId  = cal.crystal(seedId).neighbors();
      const auto nneighborsId = cal.crystal(seedId).nextNeighbors();

      double e1(seedEnergy),e9(seedEnergy),e25(seedEnergy);
      for (const auto& il : caloCrystalHitsPtrVector)
      {
          int crid = il->id();
          for (const auto& it : neighborsId)  if (it==crid) {e9 += il->energyDep(); e25 += il->energyDep(); break;}
          for (const auto& it : nneighborsId) if (it==crid) {e25 += il->energyDep(); break;}
      }

      return std::array<double,3> {e1,e9,e25};
  }



  //----------------------------------------------------------------------------------------------------
  std::array<double,4> CaloClusterFromProtoCluster::clusterTimeEnergy(const std::vector<art::Ptr<CaloCrystalHit>>& hits)
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

      return std::array<double,4>{time,timeErr,totalEnergy,totalEnergyErr};
  }

}

DEFINE_ART_MODULE(mu2e::CaloClusterFromProtoCluster);





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
