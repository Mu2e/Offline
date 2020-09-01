#include <TObject.h>
#include <TSystem.h>
#include <TFile.h>
#include "TEveEventDisplay/src/dict_classes/Collection_Filler.h"
#include "art/Framework/Principal/SubRun.h"

using namespace mu2e;
namespace mu2e{

  Collection_Filler::Collection_Filler(const Config& conf) :
    chTag_(conf.chTag()),
    crvcoinTag_(conf.crvdigiTag()),
    cluTag_(conf.cluTag()),
    cryHitTag_(conf.cryHitTag()),
    hseedTag_(conf.hseedTag()),
    kalseedTag_(conf.kalseedTag()),
    trkexttrajTag_(conf.trkexttrajTag()),
    mctrajTag_(conf.mctrajTag()),
    addHits_(conf.addHits()),
    addTracks_(conf.addTracks()),
    addClusters_(conf.addClusters()),
    addCrvHits_(conf.addCrvHits()),
    addCosmicSeedFit_(conf.addCosmicSeedFit()),
    addTrkExtTrajs_(conf.addTrkExtTrajs()),
    RecoOnly_(conf.RecoOnly()),
    FillAll_(conf.FillAll()),
    addMCTraj_(conf.addMCTraj()), 
    MCOnly_(conf.MCOnly())
  {}


  void Collection_Filler::FillRecoCollections(const art::Event& evt, Data_Collections &data, RecoDataProductName CollectionName){
    if(FillAll_ or RecoOnly_ or (addHits_ and CollectionName == ComboHits)){ 
      auto chH = evt.getValidHandle<mu2e::ComboHitCollection>(chTag_);
      data.chcol = chH.product();
    }
    if(FillAll_ or RecoOnly_ or CollectionName == CaloCrystalHits){
      auto chH = evt.getValidHandle<mu2e::CaloCrystalHitCollection>(cryHitTag_);
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
      auto chH = evt.getValidHandle<mu2e::KalSeedCollection>(kalseedTag_);
      data.kalseedcol = chH.product();
    }
    if(FillAll_ or RecoOnly_ or (addCrvHits_ and CollectionName==CRVRecoPulses)){
      auto chH = evt.getValidHandle<mu2e::CrvRecoPulseCollection>(crvcoinTag_);
      data.crvcoincol = chH.product();
    }
    if(FillAll_ or RecoOnly_ or (addTrkExtTrajs_ and CollectionName==TrkExtTrajectories)){
      auto chH = evt.getValidHandle<mu2e::TrkExtTrajCollection>(trkexttrajTag_);
      data.trkextcol = chH.product();
    }
  }


  void Collection_Filler::FillMCCollections(const art::Event& evt, Data_Collections &data, MCDataProductName CollectionName){
    if(FillAll_ or MCOnly_ or (addMCTraj_ and CollectionName == MCTrajectories)){ 
      auto chH = evt.getValidHandle<mu2e::MCTrajectoryCollection>(mctrajTag_);
      data.mctrajcol = chH.product();
    }
  }

}



