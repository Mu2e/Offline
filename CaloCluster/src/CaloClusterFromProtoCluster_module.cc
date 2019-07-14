//
// Module to produce the calorimeter clusters from proto-clusters. See MakeCaloProtocluster for the proto-cluster formation.
//
// The strategy is to attach the split-off to the main cluster first, take the closest cluster if the split-off time is compatible
// with several clusters. Then associate energetic clusters between them, including the reattached split-off of each cluster in the
// comparison.
//
// Note: the cluster center-of-gravity is calculated in the calorimeter section front face frame
//
// Original author: B. Echenard


#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "cetlib_except/exception.h"


#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CaloCluster/inc/ClusterAssociator.hh"
#include "CaloCluster/inc/ClusterMoments.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "RecoDataProducts/inc/CaloProtoClusterCollection.hh"


#include "TH1D.h"
#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <memory>
#include <tuple>
#include <array>



namespace mu2e {

  class CaloClusterFromProtoCluster : public art::EDProducer {
  public:

    explicit CaloClusterFromProtoCluster(fhicl::ParameterSet const& pset) :
      art::EDProducer{pset},
      caloClusterModuleLabel_(pset.get<std::string>("caloClusterModuleLabel")),
      mainTag_{caloClusterModuleLabel_, pset.get<std::string>("mainClusterCollName")},
      splitTag_{caloClusterModuleLabel_, pset.get<std::string>("splitClusterCollName")},
      mainToken_{consumes<CaloProtoClusterCollection>(mainTag_)},
      splitToken_{consumes<CaloProtoClusterCollection>(splitTag_)},
      deltaTime_(pset.get<double>("deltaTime")),
      maxDistSplit_(pset.get<double>("maxDistSplit")),
      maxDistMain_(pset.get<double>("maxDistMain")),
      cogTypeName_(pset.get<std::string>("cogTypeName")),
      cogType_(ClusterMoments::Linear),
      strategy_(pset.get<int>("strategy",0)),
      diagLevel_(pset.get<int>("diagLevel",0)),
      messageCategory_("CLUSTER")
    {
      produces<CaloClusterCollection>();
    }

    void beginJob() override;
    void produce(art::Event& e) override;

  private:

    std::string const caloClusterModuleLabel_;
    art::InputTag const mainTag_;
    art::InputTag const splitTag_;
    art::ProductToken<CaloProtoClusterCollection> const mainToken_;
    art::ProductToken<CaloProtoClusterCollection> const splitToken_;
    double const deltaTime_;
    double const maxDistSplit_;
    double const maxDistMain_;
    std::string const cogTypeName_;
    ClusterMoments::cogtype cogType_;
    int const strategy_;
    int const diagLevel_;
    std::string const messageCategory_;

    void makeCaloClusters(CaloClusterCollection& caloClusters,
                          const CaloProtoClusterCollection& caloClustersMain,
                          const CaloProtoClusterCollection& caloClustersSplit);

    std::array<double,3> calcEnergyLayer(const Calorimeter& cal,
                                         const std::vector<art::Ptr<CaloCrystalHit>>& caloCrystalHitsPtrVector);

    std::array<double,4> clusterTimeEnergy(const std::vector<art::Ptr<CaloCrystalHit>>& hits);
  };



  void CaloClusterFromProtoCluster::beginJob()
  {
    if (cogTypeName_.compare("Linear"))    cogType_ = ClusterMoments::Linear;
    if (cogTypeName_.compare("Logarithm")) cogType_ = ClusterMoments::Logarithm;
  }

  void CaloClusterFromProtoCluster::produce(art::Event& event)
  {
    // Check that calorimeter geometry description exists
    art::ServiceHandle<GeometryService> geom;
    if( !(geom->hasElement<Calorimeter>()) ) return;

    auto const& caloClustersMain = *event.getValidHandle(mainToken_);
    auto const& caloClustersSplit = *event.getValidHandle(splitToken_);

    //Create a new CaloCluster collection and fill it
    auto caloClusters = std::make_unique<CaloClusterCollection>();
    makeCaloClusters(*caloClusters, caloClustersMain, caloClustersSplit);

    event.put(std::move(caloClusters));
  }



  void CaloClusterFromProtoCluster::makeCaloClusters(CaloClusterCollection& caloClusters,
                                                     const CaloProtoClusterCollection& caloClustersMain,
                                                     const CaloProtoClusterCollection& caloClustersSplit)
  {
    const Calorimeter& cal = *(GeomHandle<Calorimeter>());


    ClusterAssociator associator(cal);
    associator.associateSplitOff(caloClustersMain, caloClustersSplit, deltaTime_,maxDistSplit_);



    //-- First, associate split-off to main cluster into intermediate buffer
    CaloProtoClusterCollection caloProtoClustersTemp;

    for (unsigned int imain(0); imain<caloClustersMain.size(); ++imain)
      {

        bool isSplit(false);
        std::vector<art::Ptr<CaloCrystalHit>> caloCrystalHitsPtrVector = caloClustersMain.at(imain).caloCrystalHitsPtrVector();


        //search split-offs and add their hits into the main cluster
        for (unsigned int isplit=0;isplit<caloClustersSplit.size();++isplit)
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
    for (unsigned int iproto=0;iproto<caloProtoClustersTemp.size();++iproto)
      {
        if (flagProto[iproto]) continue;

        auto        caloCrystalHitsPtrVector = caloProtoClustersTemp.at(iproto).caloCrystalHitsPtrVector();
        bool        isSplit                  = caloProtoClustersTemp.at(iproto).isSplit();
        const auto& seed                     = **caloCrystalHitsPtrVector.begin();
        int         diskId                = cal.crystal(seed.id()).diskId();

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
        std::sort(caloCrystalHitsPtrVector.begin(),caloCrystalHitsPtrVector.end(),
                  [](const art::Ptr<CaloCrystalHit>& lhs,const art::Ptr<CaloCrystalHit>& rhs) {return lhs->energyDep() > rhs->energyDep();} );


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
    std::sort(caloClusters.begin(),caloClusters.end(),[](CaloCluster& lhs, CaloCluster& rhs) {return lhs.energyDep() > rhs.energyDep();} );


  }




  //----------------------------------------------------------------------------------------------------
  std::array<double,3> CaloClusterFromProtoCluster::calcEnergyLayer(const Calorimeter& cal,
                                                                    const std::vector<art::Ptr<CaloCrystalHit>>& caloCrystalHitsPtrVector)
  {
    int seedId              = caloCrystalHitsPtrVector[0]->id();
    double seedEnergy       = caloCrystalHitsPtrVector[0]->energyDep();
    auto const neighborsId  = cal.crystal(seedId).neighbors();
    auto const nneighborsId = cal.crystal(seedId).nextNeighbors();

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
