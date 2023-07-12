// Ntuple dumper for StepPointMCs.
//
// Andrei Gaponenko, 2013

#include <string>
#include <vector>
#include <limits>
#include <cmath>

#include "cetlib_except/exception.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "TDirectory.h"
#include "TH1.h"
#include "TTree.h"

#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalSequence.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Provenance.h"
#include "art_root_io/TFileService.h"

#include "Offline/MCDataProducts/inc/StepPointMC.hh"

#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/Mu2eUtilities/inc/SimParticleGetTau.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"

#include "KinKal/General/ParticleState.hh"

namespace mu2e {

  //================================================================
  double getCharge(PDGCode::type pdgId) {
    // unlike generic conditions, MC particle data
    // should not change run-to-run, so static is safe
    // use static for efficiency
    static GlobalConstantsHandle<ParticleDataList> pdt;

    return pdt->particle(pdgId).charge();
  }

  //================================================================
  double getKineticEnergy(const StepPointMC& hit) {
    // unlike generic conditions, MC particle data
    // should not change run-to-run, so static is safe
    // use static for efficiency
    static GlobalConstantsHandle<ParticleDataList> pdt;

    const double mass = pdt->particle(hit.simParticle()->pdgId()).mass();
    return sqrt(hit.momentum().mag2() + std::pow(mass, 2)) - mass;
  }

  //================================================================
  struct VDHit {
    float x;
    float y;
    float z;
    float time;

    float px;
    float py;
    float pz;
    float pmag;
    float ek;

    float charge;
    int   pdgId;
    unsigned particleId;
    unsigned volumeCopyNumber;

    VDHit() : x(std::numeric_limits<double>::quiet_NaN())
              , y(std::numeric_limits<double>::quiet_NaN())
              , z(std::numeric_limits<double>::quiet_NaN())

              , time(std::numeric_limits<double>::quiet_NaN())

                , px(std::numeric_limits<double>::quiet_NaN())
                , py(std::numeric_limits<double>::quiet_NaN())
                , pz(std::numeric_limits<double>::quiet_NaN())
                , pmag(std::numeric_limits<double>::quiet_NaN())
                , ek(std::numeric_limits<double>::quiet_NaN())

                , charge(std::numeric_limits<double>::quiet_NaN())
                , pdgId(0)
                , particleId(-1U)
                , volumeCopyNumber(-1U)
                {}

    //----------------------------------------------------------------
    VDHit( const StepPointMC& hit)
      : x(hit.position().x())
        , y(hit.position().y())
        , z(hit.position().z())
        , time(hit.time())
        , px(hit.momentum().x())
        , py(hit.momentum().y())
        , pz(hit.momentum().z())

        , pmag(hit.momentum().mag())
          , ek(getKineticEnergy(hit))

          , charge(getCharge(hit.simParticle()->pdgId()))

          , pdgId(hit.simParticle()->pdgId())
          , particleId(hit.simParticle()->id().asUint())
          , volumeCopyNumber(hit.volumeId())
          {}

  }; // struct VDHit

  //================================================================
  class StepPointMCDumper : public art::EDAnalyzer {
    typedef std::vector<std::string> VS;
    typedef std::vector<StepPointMCCollection> VspMC;
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<std::string> hits     {Name("hitsInputTag"     ), Comment("StepPointMC collection")};
      fhicl::OptionalSequence<std::string> tauCollections     {Name("tauHitCollections"), Comment("StepPointMC collections for proper time calculation")};
      fhicl::OptionalSequence<int> decayOffCodes     {Name("decayOffPDGCodes"), Comment("decayOffPDGCodes")};
      fhicl::Atom<bool>          writeVDHit  {Name("writeVDHit"),   Comment("Write VDHit format branch"), false};
      fhicl::Atom<bool>          writeParticleState  {Name("writeParticleState"),   Comment("Write ParticleState format branch"), false};
      fhicl::Atom<bool>          writeProperTime  {Name("writeProperTime"),   Comment("Write ProperTime format branch"), false};
      fhicl::Atom<bool>          detectorSystem  {Name("detectorSystem"),   Comment("Use DetectorSystem for position information for ParticleState"), false};
    };
    typedef art::EDAnalyzer::Table<Config> Parameters;

    art::InputTag hitsInputTag_;

    bool writeVDHit_, writeParticleState_, writeProperTime_, detectorSystem_;
    VS tauHitCollections_;
    std::vector<int> decayOffCodes_;

    // Members needed to write the ntuple
    TTree *nt_;
    VDHit hit_;
    KinKal::ParticleState pstate_;
    float tau_;

    public:
    explicit StepPointMCDumper(const Parameters& pset);
    virtual void beginJob();
    virtual void analyze(const art::Event& event);
  };

  //================================================================
  StepPointMCDumper::StepPointMCDumper(const Parameters& pset)
    : art::EDAnalyzer(pset)
      , hitsInputTag_(pset().hits())
      , writeVDHit_(pset().writeVDHit())
      , writeParticleState_(pset().writeParticleState())
      , writeProperTime_(pset().writeProperTime())
      , detectorSystem_(pset().detectorSystem())
      , nt_(0)
      , tau_(0.0)
  {
    if(writeProperTime_) {
      pset().tauCollections(tauHitCollections_);
      pset().decayOffCodes(decayOffCodes_);
      // must sort to use binary_search in SimParticleGetTau
      std::sort(decayOffCodes_.begin(), decayOffCodes_.end());
    }
  }

  //================================================================
  void StepPointMCDumper::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;
    static const char branchDesc[] = "x/F:y/F:z/F:time/F:px/F:py/F:pz/F:pmag/F:ek/F:charge/F:pdgId/I:particleId/i:volumeCopy/i";
    nt_ = tfs->make<TTree>( "nt", "StepPointMCDumper ntuple");
    if(writeVDHit_)nt_->Branch("hits", &hit_, branchDesc);
    if(writeParticleState_)nt_->Branch("particle", &pstate_);
    if(writeProperTime_) { nt_->Branch("tau", &tau_, "tauNormalized/F"); }
  }

  //================================================================
  void StepPointMCDumper::analyze(const art::Event& event) {
    auto const& ptable = GlobalConstantsHandle<ParticleDataList>();
    GeomHandle<DetectorSystem> det;

    VspMC spMCColls;
    for ( const auto& iColl : tauHitCollections_ ){
      auto spColl = event.getValidHandle<StepPointMCCollection>(iColl);
      spMCColls.push_back( *spColl );
    }

    const auto& ih = event.getValidHandle<StepPointMCCollection>(hitsInputTag_);
    for(const auto& i : *ih) {

      if(writeVDHit_)hit_ = VDHit(i);
      if(writeParticleState_) {
        KinKal::VEC3 pos = detectorSystem_ ? KinKal::VEC3(det->toDetector(i.position())) : KinKal::VEC3(i.position());
        KinKal::VEC3 mom(i.momentum());
        double time = i.time();
        double mass = i.simParticle()->startMomentum().invariantMass();
        int charge = static_cast<int>(ptable->particle(i.simParticle()->pdgId()).charge());
        pstate_ = KinKal::ParticleState(pos,mom,time,mass,charge);
      }

      if(writeProperTime_) { tau_ = SimParticleGetTau::calculate(i, spMCColls, decayOffCodes_); }

      nt_->Fill();
    }

  } // analyze(event)

  //================================================================

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::StepPointMCDumper)
