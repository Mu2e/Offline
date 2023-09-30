#ifndef RobustMultiHelixFinderDiag_hh
#define RobustMultiHelixFinderDiag_hh

#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "fhiclcpp/ParameterSet.h"

#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/MCDataProducts/inc/MCRelationship.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/Mu2eUtilities/inc/ModuleHistToolBase.hh"
#include "Offline/TrkPatRec/inc/RobustMultiHelixFinder_types.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"

#include "TTree.h"

#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "art/Framework/Core/EDAnalyzer.h"
#include "Offline/GeometryService/inc/DetectorSystem.hh"



namespace mu2e {

  using RobustMultiHelixFinderTypes::Data_t;
  using RobustMultiHelixFinderTypes::Config;

  class RobustMultiHelixFinderDiag: public mu2e::ModuleHistToolBase  {

    public:

      explicit RobustMultiHelixFinderDiag(const Config& config);
      explicit RobustMultiHelixFinderDiag(const fhicl::ParameterSet& PSet);
      ~RobustMultiHelixFinderDiag() = default;

      virtual int bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) override;
      virtual int fillHistograms(void* Data, int Mode = -1) override ;


    private:

      void setupBranches();

      bool           treeInit_;
      TTree*         trkdiag_;
      Data_t*        data_;
      bool           mcdiag_;
      art::InputTag  mcdigisTag_;
  };


  RobustMultiHelixFinderDiag::RobustMultiHelixFinderDiag(const Config& config) :
    treeInit_(  false),
    mcdiag_(    config.mcDiag()),
    mcdigisTag_(config.strawDigiMCCollection())
  {}

  RobustMultiHelixFinderDiag::RobustMultiHelixFinderDiag(const fhicl::ParameterSet& pset) :
    treeInit_(  false),
    mcdiag_(    pset.get<bool>(         "MCDiag",                true)),
    mcdigisTag_(pset.get<art::InputTag>("StrawDigiMCCollection","makeSD"))
  {}



  int RobustMultiHelixFinderDiag::bookHistograms(art::ServiceHandle<art::TFileService>& Tfs)
  {
    trkdiag_ = Tfs->make<TTree>("rhfdiag","robust multi helix finder diagnostics");
    return 0;
  }

  void RobustMultiHelixFinderDiag::setupBranches()
  {
    trkdiag_->Branch("iev",       &data_->iev_,       "iev/I");
    trkdiag_->Branch("Nch",       &data_->Nch_,       "Nch/I");
    trkdiag_->Branch("chNhit",    &data_->chNhit_,    "chNhit[Nch]/I");
    trkdiag_->Branch("chPdg",     &data_->chPdg_,     "chPdg[Nch]/I");
    trkdiag_->Branch("chCrCode",  &data_->chCrCode_,  "chCrCode[Nch]/I");
    trkdiag_->Branch("chSimId",   &data_->chSimId_,   "chSimId[Nch]/I");
    trkdiag_->Branch("chTime",    &data_->chTime_,    "chTime[Nch]/F");
    trkdiag_->Branch("chPhi",     &data_->chPhi_,     "chPhi[Nch]/F");
    trkdiag_->Branch("chRad",     &data_->chRad_,     "chRad[Nch]/F");
    trkdiag_->Branch("chX",       &data_->chX_,       "chX[Nch]/F");
    trkdiag_->Branch("chY",       &data_->chY_,       "chY[Nch]/F");
    trkdiag_->Branch("chZ",       &data_->chZ_,       "chZ[Nch]/F");
    trkdiag_->Branch("chUId",     &data_->chUId_,     "chUId[Nch]/I");
    trkdiag_->Branch("chTerr",    &data_->chTerr_,    "chTerr[Nch]/F");
    trkdiag_->Branch("chWerr",    &data_->chWerr_,    "chWerr[Nch]/F");
    trkdiag_->Branch("chWDX",     &data_->chWDX_,     "chWDX[Nch]/F");
    trkdiag_->Branch("chWDY",     &data_->chWDY_,     "chWDY[Nch]/F");
    trkdiag_->Branch("chStrawId", &data_->chStrawId_);
    trkdiag_->Branch("Nhel",      &data_->Nhel_,      "Nhel/I");
    trkdiag_->Branch("helhel",    &data_->helhel_,    "helhel[Nhel]/I");
    trkdiag_->Branch("helrad",    &data_->helrad_,    "helrad[Nhel]/F");
    trkdiag_->Branch("helrcen",   &data_->helrcen_,   "helrcen[Nhel]/F");
    trkdiag_->Branch("helfcen",   &data_->helfcen_,   "helfcen[Nhel]/F");
    trkdiag_->Branch("hellam",    &data_->hellam_,    "hellam[Nhel]/F");
    trkdiag_->Branch("helchi2",   &data_->helchi2_,   "helchi2[Nhel]/F");
    trkdiag_->Branch("heldzdt",   &data_->heldzdt_,   "heldzdt[Nhel]/F");
    trkdiag_->Branch("helfz0",    &data_->helfz0_,    "helfz0[Nhel]/F");
    trkdiag_->Branch("helnhi",    &data_->helnhi_,    "helnhi[Nhel]/I");
    trkdiag_->Branch("helhits",   &data_->helhits_);
    treeInit_ = true;
  }


  int RobustMultiHelixFinderDiag::fillHistograms(void* Data, int Mode)
  {
      data_ = static_cast<Data_t*>(Data);

      // we need to set up the branches once the data_ structure has been allocated
      if (!treeInit_) setupBranches();

      if (mcdiag_){
          auto& mcdigis = *data_->event_->getValidHandle<StrawDigiMCCollection>(mcdigisTag_);

          for (int ich=0; ich <data_->Nch_ ; ++ich){
              std::vector<StrawDigiIndex> dids;
              data_->chcol_->fillStrawDigiIndices(ich,dids);
              const StrawDigiMC& mcdigi        = mcdigis.at(dids[0]);// taking 1st digi: is there a better idea??
              const art::Ptr<SimParticle>& spp = mcdigi.earlyStrawGasStep()->simParticle();
              data_->chPdg_[ich]    = spp->pdgId();
              data_->chCrCode_[ich] = spp->creationCode();
              data_->chSimId_[ich]  = spp->id().asInt();
          }
      }

      trkdiag_->Fill();
      return 0;
  }

}


DEFINE_ART_CLASS_TOOL(mu2e::RobustMultiHelixFinderDiag)
#endif
