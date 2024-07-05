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
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "art_root_io/TFileService.h"
#include "TTree.h"
#include "TH1F.h"
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
        fhicl::OptionalAtom<art::InputTag> VDStepCollection{   Name("VDStepCollection"),   Comment("Virtual Detector StepPointMC collection name") };
        fhicl::OptionalAtom<art::InputTag> IPAStepCollection{   Name("IPAStepCollection"),   Comment("IPA StepPointMC collection name") };
        fhicl::OptionalAtom<art::InputTag> STStepCollection{   Name("STStepCollection"),   Comment("StoppingTarget StepPointMC collection name") };
      };

      explicit SurfaceStepDiag(const art::EDAnalyzer::Table<Config>& config);
      virtual ~SurfaceStepDiag();
      virtual void beginJob();
      virtual void analyze(const art::Event& e);
    private:
      art::ProductToken<SurfaceStepCollection> sstag_;
      art::ProductToken<StepPointMCCollection> vdtag_, ipatag_, sttag_;
      TTree *ssdiag_;
      int evt_, subrun_, run_, sid_, sindex_;
      float edep_, path_, time_;
      XYZVectorF mom_, start_, mid_, end_;
      TH1F *nvd_, *nipa_, *nstfoil_,* nstwire_;
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
    ssdiag_->Branch("subrun",&subrun_,"subrun/I");
    ssdiag_->Branch("run",&run_,"run/I");
    ssdiag_->Branch("sid",&sid_,"sid/I");
    ssdiag_->Branch("sindex",&sindex_,"sindex/I");
    ssdiag_->Branch("edep",&edep_,"edep/F");
    ssdiag_->Branch("path",&path_,"path/F");
    ssdiag_->Branch("time",&time_,"time/F");
    ssdiag_->Branch("mom.",&mom_);
    ssdiag_->Branch("start.",&start_);
    ssdiag_->Branch("mid.",&mid_);
    ssdiag_->Branch("end.",&end_);
    nvd_ = tfs->make<TH1F>("nvd","N VD SurfaceSteps",100,-0.5,99.5);
    nipa_ = tfs->make<TH1F>("nipa","N IPA SurfaceSteps",100,-0.5,99.5);
    nstfoil_ = tfs->make<TH1F>("nstfoil","N ST Foil SurfaceSteps",100,-0.5,99.5);
    nstwire_ = tfs->make<TH1F>("nstwire","N ST Wire SurfaceSteps",100,-0.5,99.5);
  }

  void SurfaceStepDiag::analyze(const art::Event& evt ) {
    auto sscolH = evt.getValidHandle(sstag_);
    auto const& sscol = *sscolH.product();
    evt_ = evt.id().event();
    subrun_ = evt.id().subRun();
    run_ = evt.id().run();
    unsigned nvd(0), nipa(0), nstfoil(0), nstwire(0);
    for(auto const& ss : sscol) {
      sid_ = ss.surfaceId().id();
      sindex_ = ss.surfaceId().index();
      edep_ = ss.energyDeposit();
      path_ = ss.pathLength();
      time_ = ss.time();
      mom_ = ss.momentum();
      start_ = ss.startPosition();
      mid_ = ss.midPosition();
      end_ = ss.endPosition();
      ssdiag_->Fill();
      if(ss.surfaceId().id() < SurfaceIdDetail::IPA)nvd++;
      if(ss.surfaceId().id() == SurfaceIdDetail::IPA)nipa++;
      if(ss.surfaceId().id() == SurfaceIdDetail::ST_Foils)nstfoil++;
      if(ss.surfaceId().id() == SurfaceIdDetail::ST_Wires)nstwire++;
    }
    nvd_->Fill(nvd);
    nipa_->Fill(nipa);
    nstfoil_->Fill(nstfoil);
    nstwire_->Fill(nstwire);
  }
}
// Part of the magic that makes this class a module.
using mu2e::SurfaceStepDiag;
DEFINE_ART_MODULE(SurfaceStepDiag)
