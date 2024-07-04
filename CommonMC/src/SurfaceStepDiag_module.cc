//
// SurfaceStep diags
//
// Original author D. Brown
//
// framework
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "Offline/DataProducts/inc/GenVector.hh"
#include "art_root_io/TFileService.h"
#include "TTree.h"
// data
#include "Offline/MCDataProducts/inc/SurfaceStep.hh"

namespace mu2e
{
  class SurfaceStepDiag : public art::EDAnalyzer {
    public:
      using SPM = std::map<art::Ptr<SimParticle>,unsigned>;

      struct Config {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;
        fhicl::Atom<art::InputTag> SurfaceStepCollection{   Name("SurfaceStepCollection"),   Comment("SurfaceStep collection name") };
      };

      explicit SurfaceStepDiag(const art::EDAnalyzer::Table<Config>& config);
      virtual ~SurfaceStepDiag();
      virtual void beginJob();
      virtual void analyze(const art::Event& e);
    private:
      art::ProductToken<SurfaceStepCollection> sstag_;
      TTree *ssdiag_;
      int evt_, sid_, sindex_;
      float edep_, path_;
      XYZVectorF mom_, pos_;
  };

  SurfaceStepDiag::SurfaceStepDiag(const art::EDAnalyzer::Table<Config>& config) :
    art::EDAnalyzer(config),
    sstag_{ consumes<SurfaceStepCollection>(config().SurfaceStepCollection() ) }
  {}

  SurfaceStepDiag::~SurfaceStepDiag(){}

  void SurfaceStepDiag::beginJob() {
    // create diagnostics if requested
    art::ServiceHandle<art::TFileService> tfs;
    // detailed diagnostics
    ssdiag_=tfs->make<TTree>("ssdiag","SurfaceStep diagnostics");
    ssdiag_->Branch("evt",&evt_,"evt/I");
    ssdiag_->Branch("sid",&sid_,"sid/I");
    ssdiag_->Branch("sindex",&sindex_,"sindex/I");
    ssdiag_->Branch("edep",&edep_,"edep/F");
    ssdiag_->Branch("path",&path_,"path/F");
    ssdiag_->Branch("mom.",&mom_);
    ssdiag_->Branch("pos.",&pos_);
  }

  void SurfaceStepDiag::analyze(const art::Event& evt ) {
    auto sscolH = evt.getValidHandle(sstag_);
    auto const& sscol = *sscolH.product();
    evt_ = evt.id().event();  // add event id
    for(auto const& ss : sscol) {
      sid_ = ss.surfaceId().id();
      sindex_ = ss.surfaceId().index();
      edep_ = ss.energyDeposit();
      path_ = ss.pathLength();
      mom_ = ss.momentum();
      pos_ = ss.startPosition();
      ssdiag_->Fill();
    }
  }
}
// Part of the magic that makes this class a module.
using mu2e::SurfaceStepDiag;
DEFINE_ART_MODULE(SurfaceStepDiag)
