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
#include "art_root_io/TFileService.h"
#include "TTree.h"
#include "TH1F.h"
// data
#include "Offline/MCDataProducts/inc/SurfaceStep.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"

namespace mu2e
{
  class SurfaceStepDiag : public art::EDAnalyzer {
    public:
      using SPM = std::map<art::Ptr<SimParticle>,unsigned>;

      struct Config {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;
        fhicl::Atom<art::InputTag> SurfaceStepCollection{   Name("SurfaceStepCollection"),   Comment("SurfaceStep collection name") };
        fhicl::Atom<art::InputTag> vdstepmcs { Name("VDStepPointMCs"), Comment("Virtual Detector StepPointMC collection")};
        fhicl::Atom<art::InputTag> ipastepmcs { Name("IPAStepPointMCs"), Comment("IPAStepPointMC collection")};
        fhicl::Atom<art::InputTag> ststepmcs { Name("TargetStepPointMCs"), Comment("Stopping target StepPointMC collection")};
      };

      explicit SurfaceStepDiag(const art::EDAnalyzer::Table<Config>& config);
      virtual ~SurfaceStepDiag();
      virtual void beginJob();
      virtual void analyze(const art::Event& e);
    private:
      art::ProductToken<SurfaceStepCollection> surfsteps_;
      art::ProductToken<StepPointMCCollection> vdstepmcs_, ipastepmcs_, ststepmcs_;
      TTree *ssdiag_;
      int evt_, subrun_, run_, sid_, sindex_;
      float edep_, path_, time_;
      XYZVectorF mom_, start_, mid_, end_;
      TH1F *nvdspmc_, *nipaspmc_, *nstspmc_;
      TH1F *nvd_, *nipa_, *nstfoil_,* nstwire_;
  };

  SurfaceStepDiag::SurfaceStepDiag(const art::EDAnalyzer::Table<Config>& config) :
    art::EDAnalyzer(config),
    surfsteps_{ consumes<SurfaceStepCollection>(config().SurfaceStepCollection() ) },
    vdstepmcs_{ consumes<StepPointMCCollection>(config().vdstepmcs())},
    ipastepmcs_{ consumes<StepPointMCCollection>(config().ipastepmcs())},
    ststepmcs_{ consumes<StepPointMCCollection>(config().ststepmcs())}
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
    nvdspmc_ = tfs->make<TH1F>("nvdspmc","N VD StepPointMCs",100,-0.5,99.5);
    nipaspmc_ = tfs->make<TH1F>("nipaspmc","N IPA StepPointMCs",100,-0.5,99.5);
    nstspmc_ = tfs->make<TH1F>("nstspmc","N ST StepPointMCs",100,-0.5,99.5);
    nvd_ = tfs->make<TH1F>("nvd","N VD SurfaceSteps",100,-0.5,99.5);
    nipa_ = tfs->make<TH1F>("nipa","N IPA SurfaceSteps",100,-0.5,99.5);
    nstfoil_ = tfs->make<TH1F>("nstfoil","N ST Foil SurfaceSteps",100,-0.5,99.5);
    nstwire_ = tfs->make<TH1F>("nstwire","N ST Wire SurfaceSteps",100,-0.5,99.5);
  }

  void SurfaceStepDiag::analyze(const art::Event& event ) {
// first histogram StepPointMC precursor sizes
    auto const& vdspmccol_h = event.getValidHandle<StepPointMCCollection>(vdstepmcs_);
    auto const& vdspmccol = *vdspmccol_h;
    nvdspmc_->Fill(vdspmccol.size());
    auto const& ipaspmccol_h = event.getValidHandle<StepPointMCCollection>(ipastepmcs_);
    auto const& ipaspmccol = *ipaspmccol_h;
    nipaspmc_->Fill(ipaspmccol.size());
    auto const& stspmccol_h = event.getValidHandle<StepPointMCCollection>(ststepmcs_);
    auto const& stspmccol = *stspmccol_h;
    nstspmc_->Fill(stspmccol.size());


    auto sscolH = event.getValidHandle(surfsteps_);
    auto const& sscol = *sscolH.product();
    evt_ = event.id().event();
    subrun_ = event.id().subRun();
    run_ = event.id().run();
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
