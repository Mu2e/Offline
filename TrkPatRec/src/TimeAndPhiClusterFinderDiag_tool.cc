#ifndef TimeAndPhiClusterFinderDiag_hh
#define TimeAndPhiClusterFinderDiag_hh

#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "fhiclcpp/ParameterSet.h"

#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/Mu2eUtilities/inc/ModuleHistToolBase.hh"
#include "Offline/TrkPatRec/inc/TimeAndPhiClusterFinder_types.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"

#include "TTree.h"


#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "art/Framework/Core/EDAnalyzer.h"
#include "Offline/GeometryService/inc/DetectorSystem.hh"



namespace mu2e {

  using TimeAndPhiClusterFinderTypes::Data_t;
  using TimeAndPhiClusterFinderTypes::Config;

  class TimeAndPhiClusterFinderDiag: public mu2e::ModuleHistToolBase  {

    public:

      explicit TimeAndPhiClusterFinderDiag(const Config& config);
      explicit TimeAndPhiClusterFinderDiag(const fhicl::ParameterSet& PSet);
      ~TimeAndPhiClusterFinderDiag() = default;

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


  TimeAndPhiClusterFinderDiag::TimeAndPhiClusterFinderDiag(const Config& config) :
    treeInit_(  false),
    mcdiag_(    config.mcDiag()),
    mcdigisTag_(config.strawDigiMCCollection())
  {}

  TimeAndPhiClusterFinderDiag::TimeAndPhiClusterFinderDiag(const fhicl::ParameterSet& pset) :
    treeInit_(  false),
    mcdiag_(    pset.get<bool>(         "MCDiag",                true)),
    mcdigisTag_(pset.get<art::InputTag>("StrawDigiMCCollection","makeSD"))
  {}



  int TimeAndPhiClusterFinderDiag::bookHistograms(art::ServiceHandle<art::TFileService>& Tfs)
  {
    trkdiag_ = Tfs->make<TTree>("tpcdiag","time and phi cluster diagnostics");
    return 0;
  }

  void TimeAndPhiClusterFinderDiag::setupBranches()
  {
    trkdiag_->Branch("iev",       &data_->iev_,       "iev/I");
    trkdiag_->Branch("Nch",       &data_->Nch_,       "Nch/I");
    trkdiag_->Branch("chSel",     &data_->chSel_,     "chSel[Nch]/I");
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
    trkdiag_->Branch("nhit1",     &data_->nhit1_,     "nhit1/I");
    trkdiag_->Branch("hitIdx1",   &data_->hitIdx1_,   "hitIdx1[nhit1]/I");
    trkdiag_->Branch("nclu1",     &data_->nclu1_,     "nclu1[nhit1]/I");
    trkdiag_->Branch("nhit2",     &data_->nhit2_,     "nhit2/I");
    trkdiag_->Branch("hitIdx2",   &data_->hitIdx2_,   "hitIdx2[nhit2]/I");
    trkdiag_->Branch("nclu2",     &data_->nclu2_,     "nclu2[nhit2]/I");

    treeInit_ = true;
  }


  int TimeAndPhiClusterFinderDiag::fillHistograms(void* Data, int Mode)
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


DEFINE_ART_CLASS_TOOL(mu2e::TimeAndPhiClusterFinderDiag)
#endif
