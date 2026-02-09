#include "Offline/TEveEventDisplay/src/dict_classes/Collection_Filler.h"
using namespace mu2e;
namespace mu2e{

#ifndef __ROOTCLING__
  Collection_Filler::Collection_Filler(const Config& conf) :
    chTag_(conf.chTag()),
    tcTag_(conf.tcTag()),
    crvcoinTag_(conf.crvdigiTag()),
    cluTag_(conf.cluTag()),
    cryHitTag_(conf.cryHitTag()),
    hseedTag_(conf.hseedTag()),
    kalseedTag_(conf.kalseedTag()),
    mctrajTag_(conf.mctrajTag()),
    addHits_(conf.addHits()),
    addTimeClusters_(conf.addTimeClusters()),
    addTrkHits_(conf.addTrkHits()),
    addTracks_(conf.addTracks()),
    addClusters_(conf.addClusters()),
    addCrvHits_(conf.addCrvHits()),
    addCosmicSeedFit_(conf.addCosmicSeedFit()),
    RecoOnly_(conf.RecoOnly()),
    FillAll_(conf.FillAll()),
    addMCTraj_(conf.addMCTraj()),
    addKKTracks_(conf.addKKTracks()),
    MCOnly_(conf.MCOnly())
  {}
#endif

  /*------------Function to turn InputTag to string for track labels:-------------*/
  template <typename T>
  std::string TurnNameToString( const T& value )
  {
    std::ostringstream ss;
    ss << value;
    return ss.str();
  }


  /*------------Function to fill RecoDataProduct lists:-------------*/
  void Collection_Filler::FillRecoCollections(const art::Event& evt, Data_Collections &data, RecoDataProductName CollectionName){
    if(FillAll_ or RecoOnly_ or (addHits_ and CollectionName == ComboHits)){
      auto chH = evt.getValidHandle<mu2e::ComboHitCollection>(chTag_);
      data.chcol = chH.product();
    }
    if(FillAll_ or RecoOnly_ or (addTimeClusters_ and CollectionName == TimeClusters)){
      auto chH = evt.getValidHandle<mu2e::TimeClusterCollection>(tcTag_);
      data.tccol = chH.product();
    }
    if(FillAll_ or RecoOnly_ or (addTrkHits_ and CollectionName == ComboHits)){
      auto chH = evt.getValidHandle<mu2e::ComboHitCollection>(chTag_);
      data.chcol = chH.product();
    }
    if(FillAll_ or RecoOnly_ or CollectionName == CaloCrystalHits){
      auto chH = evt.getValidHandle<mu2e::CaloHitCollection>(cryHitTag_);
      data.cryHitcol = chH.product();
    }
    if(FillAll_ or RecoOnly_ or (addClusters_ and CollectionName==CaloClusters)){
      auto chH = evt.getValidHandle<mu2e::CaloClusterCollection>(cluTag_);
      data.clustercol = chH.product();
    }
    if(FillAll_ or RecoOnly_ or (addCosmicSeedFit_ and CollectionName==CosmicTracks)){
      auto chH = evt.getValidHandle<mu2e::CosmicTrackSeedCollection>("CosmicTrackFinderTimeFit");
      data.cosmiccol = chH.product();
    }
    if(FillAll_ or RecoOnly_ or (addTracks_ and CollectionName==HelixSeeds)){
      auto chH = evt.getValidHandle<mu2e::HelixSeedCollection>(hseedTag_);
      data.hseedcol = chH.product();
    }
    if(FillAll_ or RecoOnly_ or (addTracks_ and CollectionName==KalSeeds)){

          for(const auto &tag : kalseedTag_){
            auto chH = evt.getValidHandle<mu2e::KalSeedCollection>(tag);
            data.kalseedcol = chH.product();
            data.track_list.push_back(data.kalseedcol);

            std::string name = TurnNameToString(tag);
            std::cout<<"Plotting KalSeed Instance: "<<name<<std::endl;
            data.track_labels.push_back(name);

          }
          data.track_tuple = std::make_tuple(data.track_labels,data.track_list);
    }
    if(FillAll_ or RecoOnly_ or (addCrvHits_ and CollectionName==CRVRecoPulses)){
      auto chH = evt.getValidHandle<mu2e::CrvRecoPulseCollection>(crvcoinTag_);
      data.crvcoincol = chH.product();
    }
  }

  /*------------Function to fill MCDataProduct lists:-------------*/
  void Collection_Filler::FillMCCollections(const art::Event& evt, Data_Collections &data, MCDataProductName CollectionName){
    if(FillAll_ or MCOnly_ or (addMCTraj_ and CollectionName == MCTrajectories)){
      auto chH = evt.getValidHandle<mu2e::MCTrajectoryCollection>(mctrajTag_);
      std::string name = TurnNameToString(mctrajTag_);
      std::cout<<"Plotting MCtraj Instance: "<<name<<std::endl;
      data.mctrajcol = chH.product();
    }
  }

}
