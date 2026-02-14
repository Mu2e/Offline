#ifndef Collection_Filler_h
#define Collection_Filler_h
//ROOT
#include <TObject.h>
#include <TROOT.h>
#include <TGComboBox.h>
#include <TGListBox.h>
#include <iostream>
#include <vector>
#include <TSystem.h>
#include <TFile.h>
//Cosmics:
#include "Offline/RecoDataProducts/inc/CosmicTrackSeed.hh"
//Calo:
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
//MC Products:
#include "Offline/MCDataProducts/inc/MCTrajectoryCollection.hh"
//Kalman Tracks
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
//Tracker Hits:
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
//CRV:
#include "Offline/RecoDataProducts/inc/CrvRecoPulse.hh"
//#include "RecoDataProducts/inc/CrvCoincidenceCluster.hh"
//Art:
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#ifndef __ROOTCLING__
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#endif

using namespace CLHEP;

#include "Offline/TEveEventDisplay/src/dict_classes/Data_Collections.h"

namespace mu2e{

  enum RecoDataProductName {ComboHits, TimeClusters, CaloCrystalHits, CaloClusters, CosmicTracks, HelixSeeds, KalSeeds, CRVRecoPulses};
  enum MCDataProductName {MCTrajectories};

        class Collection_Filler
        {
  public:
#ifndef __ROOTCLING__

    struct Config{
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      //RecoData Products
      fhicl::Atom<int> diagLevel{Name("diagLevel"), Comment("for info"),0};
      fhicl::Atom<art::InputTag>chTag{Name("ComboHitCollection"),Comment("chTag"), "makePH"};
      fhicl::Atom<art::InputTag>tcTag{Name("TimeClusterCollection"),Comment("tcTag"), "makePH"};
      fhicl::Atom<art::InputTag>crvdigiTag{Name("CrvRecoPulseCollection"),Comment("crvTag")};
      fhicl::Atom<art::InputTag>cosmicTag{Name("CosmicTrackSeedCollection"),Comment("cosmicTag"), "CosmicTrackFinderTimeFit"};
      fhicl::Atom<art::InputTag>cluTag{Name("CaloClusterCollection"),Comment("cluTag")};
      fhicl::Atom<art::InputTag>cryHitTag{Name("CaloHitCollection"),Comment("cryHitTag")};
      fhicl::Atom<art::InputTag>hseedTag{Name("HelixSeedCollection"),Comment("hseedTag")};
      fhicl::Sequence<art::InputTag>kalseedTag{Name("KalSeedCollection"),Comment("kalseedTag")};
      //MC Data Products
      fhicl::Atom<art::InputTag>mctrajTag{Name("MCTrajectoryCollection"),Comment("mctrajTag")};

      //To add RecoDataProducts
      fhicl::Atom<bool> addHits{Name("addHits"), Comment("set to add the hits"),false};
      fhicl::Atom<bool> addTimeClusters{Name("addTimeClusters"), Comment("set to add the TC hits"),false};
      fhicl::Atom<bool> addTrkHits{Name("addTrkHits"), Comment("set to add the Trk hits"),false};
      fhicl::Atom<bool> addTracks{Name("addTracks"), Comment("set to add tracks"),false};
      fhicl::Atom<bool> addClusters{Name("addClusters"), Comment("set to add calo lusters"),false};
      fhicl::Atom<bool> addCrvHits{Name("addCrvHits"), Comment("set to add crv hits"),false};
      fhicl::Atom<bool> addCrystallHits{Name("addCrystalHits"), Comment("for calo cry hits"), false};
      fhicl::Atom<bool> addCosmicSeedFit{Name("addCosmicSeedFit"), Comment("for fitted cosmic track"), false};
      fhicl::Atom<bool> RecoOnly{Name("RecoOnly"), Comment("set to see only Reco Data Products"), false};
      fhicl::Atom<bool> FillAll{Name("FillAll"), Comment("to see all available products"), false};
      fhicl::Atom<bool> addMCTraj{Name("addMCTraj"), Comment("set to add MC trajectories"), false};
      fhicl::Atom<bool> addKKTracks{Name("addKKTracks"), Comment("set to add KinKal traj"), false};
      fhicl::Atom<bool> MCOnly{Name("MCOnly"), Comment("set to see only MC Data Products"), false};
    };

    explicit Collection_Filler(const Config& conf);
#endif
    Collection_Filler(const Collection_Filler &);
    Collection_Filler& operator=(const Collection_Filler &);

    //RecoDataProducts:
    art::InputTag chTag_;
    art::InputTag tcTag_;
    art::InputTag crvcoinTag_;
    art::InputTag cosmicTag_;
    art::InputTag cluTag_;
    art::InputTag cryHitTag_;
    art::InputTag hseedTag_;
    std::vector<art::InputTag> kalseedTag_;
    art::InputTag trkexttrajTag_;

    //MCDataProdutcs:
    art::InputTag mctrajTag_;

    art::Event *_event;
    art::Run *_run;
    bool addHits_, addTimeClusters_, addTrkHits_, addTracks_, addClusters_, addCrvHits_, addCosmicSeedFit_, isCosmic_, RecoOnly_,  FillAll_, addMCCaloDigis_, addMCTraj_, addKKTracks_, MCOnly_;
    void FillRecoCollections(const art::Event& evt, Data_Collections &data, RecoDataProductName code);
    void FillMCCollections(const art::Event& evt, Data_Collections &data, MCDataProductName code);
    virtual ~Collection_Filler(){};

  private:
#ifndef __ROOTCLING__
    Config _conf;
    #endif
  ClassDef(Collection_Filler,0);
};

}

#endif
