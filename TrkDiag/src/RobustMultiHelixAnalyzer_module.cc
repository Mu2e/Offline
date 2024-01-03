#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "canvas/Utilities/InputTag.h"
#include "art_root_io/TFileService.h"

#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/MCDataProducts/inc/MCRelationship.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/RecoDataProducts/inc/RobustHelix.hh"
#include "Offline/RecoDataProducts/inc/StrawHitPosition.hh"

#include "TTree.h"


namespace mu2e {

  class RobustMultiHelixAnalyzer : public art::EDAnalyzer
  {
    public:
      struct Config {
        using Name    = fhicl::Name;
        using Comment = fhicl::Comment;
        fhicl::Atom<art::InputTag> strawDigiMCCollection    {Name("StrawDigiMCCollection"), Comment("StrawDigiMC collection name")                 };
        fhicl::Atom<art::InputTag> comboHitCollection       {Name("ComboHitCollection"),    Comment("ComboHit collection name")                    };
        fhicl::Atom<art::InputTag> negHelixSeedCollection   {Name("NegHelixSeedCollection"),Comment("negative helicity HelixSeed collection name") };
        fhicl::Atom<art::InputTag> posHelixSeedCollection   {Name("PosHelixSeedCollection"),Comment("positive helicity HelixSeed collection name") };
        fhicl::Atom<bool>          mcDiag                   {Name("MCDiag"),                Comment("Switch to enable MC diag")                    };
      };
      explicit RobustMultiHelixAnalyzer(const art::EDAnalyzer::Table<Config>& conf);

      void beginJob() override;
      void analyze(const art::Event& evt) override;


    private:
      art::InputTag  mcdigisTag_;
      art::InputTag  comboHitsTag_;
      art::InputTag  negHelixSeedTag_;
      art::InputTag  posHelixSeedTag_;
      bool           mcdiag_;
      TTree*         trkdiag_;

      enum  {kMaxHits = 8192, kMaxHelix = 128};
      Int_t   iev_;
      Int_t   Nch_,chNhit_[kMaxHits],chPdg_[kMaxHits],chCrCode_[kMaxHits],chSimId_[kMaxHits],chUId_[kMaxHits];
      Float_t chTime_[kMaxHits],chPhi_[kMaxHits],chRad_[kMaxHits],chX_[kMaxHits],chY_[kMaxHits],chZ_[kMaxHits];
      Float_t chTerr_[kMaxHits],chWerr_[kMaxHits],chWDX_[kMaxHits],chWDY_[kMaxHits];

      Int_t   Nhel_,helhel_[kMaxHelix],helnhi_[kMaxHelix];
      Float_t helrad_[kMaxHelix],helrcen_[kMaxHelix],helfcen_[kMaxHelix],hellam_[kMaxHelix],helfz0_[kMaxHelix];
      Float_t helchi2_[kMaxHelix],heldzdt_[kMaxHelix];
      std::vector<std::vector<int>>      helhits_;
      std::vector<std::vector<uint16_t>> chStrawId_;
  };





  //================================================================
  RobustMultiHelixAnalyzer::RobustMultiHelixAnalyzer(const art::EDAnalyzer::Table<Config>& config) :
    art::EDAnalyzer(config),
    mcdigisTag_     (config().strawDigiMCCollection()),
    comboHitsTag_   (config().comboHitCollection()),
    negHelixSeedTag_(config().negHelixSeedCollection()),
    posHelixSeedTag_(config().posHelixSeedCollection()),
    mcdiag_         (config().mcDiag()),
    trkdiag_        ()
  {}

  //================================================================
  void RobustMultiHelixAnalyzer::beginJob()
  {
    art::ServiceHandle<art::TFileService> tfs;
    trkdiag_ = tfs->make<TTree>("rhfdiag","robust multi helix finder diagnostics");
    trkdiag_->Branch("iev",       &iev_,       "iev/I");
    trkdiag_->Branch("Nch",       &Nch_,       "Nch/I");
    trkdiag_->Branch("chNhit",    &chNhit_,    "chNhit[Nch]/I");
    trkdiag_->Branch("chPdg",     &chPdg_,     "chPdg[Nch]/I");
    trkdiag_->Branch("chCrCode",  &chCrCode_,  "chCrCode[Nch]/I");
    trkdiag_->Branch("chSimId",   &chSimId_,   "chSimId[Nch]/I");
    trkdiag_->Branch("chTime",    &chTime_,    "chTime[Nch]/F");
    trkdiag_->Branch("chPhi",     &chPhi_,     "chPhi[Nch]/F");
    trkdiag_->Branch("chRad",     &chRad_,     "chRad[Nch]/F");
    trkdiag_->Branch("chX",       &chX_,       "chX[Nch]/F");
    trkdiag_->Branch("chY",       &chY_,       "chY[Nch]/F");
    trkdiag_->Branch("chZ",       &chZ_,       "chZ[Nch]/F");
    trkdiag_->Branch("chUId",     &chUId_,     "chUId[Nch]/I");
    trkdiag_->Branch("chTerr",    &chTerr_,    "chTerr[Nch]/F");
    trkdiag_->Branch("chWerr",    &chWerr_,    "chWerr[Nch]/F");
    trkdiag_->Branch("chWDX",     &chWDX_,     "chWDX[Nch]/F");
    trkdiag_->Branch("chWDY",     &chWDY_,     "chWDY[Nch]/F");
    trkdiag_->Branch("chStrawId", &chStrawId_);
    trkdiag_->Branch("Nhel",      &Nhel_,      "Nhel/I");
    trkdiag_->Branch("helhel",    &helhel_,    "helhel[Nhel]/I");
    trkdiag_->Branch("helrad",    &helrad_,    "helrad[Nhel]/F");
    trkdiag_->Branch("helrcen",   &helrcen_,   "helrcen[Nhel]/F");
    trkdiag_->Branch("helfcen",   &helfcen_,   "helfcen[Nhel]/F");
    trkdiag_->Branch("hellam",    &hellam_,    "hellam[Nhel]/F");
    trkdiag_->Branch("helchi2",   &helchi2_,   "helchi2[Nhel]/F");
    trkdiag_->Branch("heldzdt",   &heldzdt_,   "heldzdt[Nhel]/F");
    trkdiag_->Branch("helfz0",    &helfz0_,    "helfz0[Nhel]/F");
    trkdiag_->Branch("helnhi",    &helnhi_,    "helnhi[Nhel]/I");
    trkdiag_->Branch("helhits",   &helhits_);
  }

  //================================================================
  void RobustMultiHelixAnalyzer::analyze(const art::Event& event) {

      auto chH = event.getValidHandle<ComboHitCollection>(comboHitsTag_);
      const ComboHitCollection& chcol(*chH);

      auto nhH = event.getValidHandle<HelixSeedCollection>(negHelixSeedTag_);
      const HelixSeedCollection& negHelices(*nhH);

      auto phH = event.getValidHandle<HelixSeedCollection>(posHelixSeedTag_);
      const HelixSeedCollection& posHelices(*phH);

      Nch_=Nhel_=0;
      chStrawId_.clear();
      helhits_.clear();

      iev_ = event.id().event();

      Nch_ = chcol.size();
      for (unsigned ich=0; ich<chcol.size();++ich)
      {
         chTime_[ich] = chcol[ich].correctedTime();
         chX_[ich]    = chcol[ich].pos().x();
         chY_[ich]    = chcol[ich].pos().y();
         chZ_[ich]    = chcol[ich].pos().z();
         chRad_[ich]  = sqrt(chX_[ich]*chX_[ich]+chY_[ich]*chY_[ich]);
         chPhi_[ich]  = chcol[ich].pos().phi();
         chNhit_[ich] = chcol[ich].nStrawHits();
         chUId_[ich]  = chcol[ich].strawId().uniquePanel();
         chTerr_[ich] = chcol[ich].posRes(StrawHitPosition::trans);
         chWerr_[ich] = chcol[ich].posRes(StrawHitPosition::wire);
         chWDX_[ich]  = chcol[ich].uDir2D().x();
         chWDY_[ich]  = chcol[ich].uDir2D().y();

         std::vector<StrawHitIndex> strawHitIdxs;
         chcol.fillStrawHitIndices(ich,strawHitIdxs);
         chStrawId_.push_back(strawHitIdxs);
      }

      if (mcdiag_){
        auto& mcdigis = *event.getValidHandle<StrawDigiMCCollection>(mcdigisTag_);
        for (size_t ich=0; ich<chcol.size(); ++ich){
          std::vector<StrawDigiIndex> dids;
          chcol.fillStrawDigiIndices(ich,dids);
          const StrawDigiMC& mcdigi        = mcdigis.at(dids[0]);
          const art::Ptr<SimParticle>& spp = mcdigi.earlyStrawGasStep()->simParticle();
          chPdg_[ich]    = spp->pdgId();
          chCrCode_[ich] = spp->creationCode();
          chSimId_[ich]  = spp->id().asInt();
        }
      }

      for (const auto& helix : posHelices){
        std::vector<int> hhits;
        for (auto hi : helix.hits()){
          for (size_t j=0;j<chcol.size();++j){
            if (hi.index(0)==chcol[j].index(0)) {hhits.push_back(j); break;}
          }
        }
        helhel_[Nhel_]  = 1;
        helrad_[Nhel_]  = helix.helix().radius();
        helrcen_[Nhel_] = helix.helix().rcent();
        helfcen_[Nhel_] = helix.helix().fcent();
        hellam_[Nhel_]  = helix.helix().lambda();
        helfz0_[Nhel_]  = helix.helix().fz0();
        helchi2_[Nhel_] = helix.helix().chi2dZPhi();
        heldzdt_[Nhel_] = helix.helix().chi2dXY();
        helnhi_[Nhel_]  = hhits.size();
        helhits_.push_back(hhits);
        ++Nhel_;
     }

     for (const auto& helix : negHelices){
        std::vector<int> hhits;
        for (auto hi : helix.hits()){
          for (size_t j=0;j<chcol.size();++j){
            if (hi.index(0)==chcol[j].index(0)) {hhits.push_back(j); break;}
          }
        }
        helhel_[Nhel_]  = -1;
        helrad_[Nhel_]  = helix.helix().radius();
        helrcen_[Nhel_] = helix.helix().rcent();
        helfcen_[Nhel_] = helix.helix().fcent();
        hellam_[Nhel_]  = helix.helix().lambda();
        helfz0_[Nhel_]  = helix.helix().fz0();
        heldzdt_[Nhel_] = helix.helix().chi2dXY();
        helnhi_[Nhel_]  = hhits.size();
        helchi2_[Nhel_] = helix.helix().chi2dZPhi();
        helhits_.push_back(hhits);
        ++Nhel_;
      }

      trkdiag_->Fill();
  }
}
DEFINE_ART_MODULE(mu2e::RobustMultiHelixAnalyzer)
