// Write into an ntuple information about time, position, and
// (optionally) proper time of SimParticle end points.  There are two
// ways to specify the set of particles to process:
//
// 1) dumpSimParticleLeaves=false (default), inputCollection is a
//    SimParticlePtrCollection that explicitly lists what to dump.
//
// 2) dumpSimParticleLeaves=true, inputCollection is a SimParticle
//    collection.  The leaves of the SimParticle tree will be dumped.
//
// Andrei Gaponenko, 2013, 2015

#include <string>
#include <algorithm>

#include "cetlib_except/exception.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/SimParticlePtrCollection.hh"
#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "Mu2eUtilities/inc/SimParticleGetTau.hh"

#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/PhysicsParams.hh"

#include "TTree.h"

#include <algorithm>
#include <iterator>

namespace mu2e {

  namespace {

    typedef std::vector<std::string> VS;
    typedef std::vector<StepPointMCCollection> VspMC;

    // This should be minimal info, we'll want to load this
    // in memory in consumer jobs.  This is NOT an analysis ntuple!
    struct StopInfo {
      float x;
      float y;
      float z;
      float t;
      float tau; // proper time, for stopped pion weights

      StopInfo() : x(), y(), z(), t(), tau() {}

      StopInfo(const art::Ptr<SimParticle>& p, const VspMC& spMCcolls, float tt)
        : x(p->endPosition().x())
        , y(p->endPosition().y())
        , z(p->endPosition().z())
        , t(p->endGlobalTime())
        , tau(tt)
      {
        if(!p->endDefined()) {
          throw cet::exception("BADINPUTS")
            <<"StoppedParticlesDumper: input SimParticle does not have end defined!\n";
        }
      }
    };

  }// namespace

  //================================================================
  class StoppedParticlesDumper : public art::EDAnalyzer {
  public:
    explicit StoppedParticlesDumper(fhicl::ParameterSet const& pset);
    void beginJob() override;
    void analyze(const art::Event& evt) override;
  private:
    bool dumpSimParticleLeaves_;
    art::InputTag input_;
    bool writeProperTime_;
    VS hitColls_;

    std::vector<int> decayOffCodes_;

    TTree *nt_;
    StopInfo data_;

    bool is_leave(const SimParticle& p);
    void process(const art::Ptr<SimParticle>& p, const VspMC& spMCcolls);
  };

  //================================================================
  StoppedParticlesDumper::StoppedParticlesDumper(const fhicl::ParameterSet& pset) :
    art::EDAnalyzer(pset),
    dumpSimParticleLeaves_(pset.get<bool>("dumpSimParticleLeaves", false)),
    input_(pset.get<std::string>("inputCollection")),
    writeProperTime_(pset.get<bool>("writeProperTime")),
    hitColls_( writeProperTime_ ? pset.get<VS>("hitCollections") : VS() ),
    nt_()
  {
    if(writeProperTime_) {
      decayOffCodes_ = pset.get<std::vector<int> >("decayOffPDGCodes");
      // must sort to use binary_search in SimParticleGetTau
      std::sort(decayOffCodes_.begin(), decayOffCodes_.end());
    }
  }

  //================================================================
  void StoppedParticlesDumper::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;
    std::string branchDesc("x/F:y/F:z/F:time/F");
    if(writeProperTime_) {
      branchDesc += ":tauNormalized/F";
    }
    nt_ = tfs->make<TTree>( "stops", "Stopped particles ntuple");
    nt_->Branch("stops", &data_, branchDesc.c_str());
  }

  //================================================================
  void StoppedParticlesDumper::analyze(const art::Event& event) {

    VspMC spMCColls;

    for ( const auto& iColl : hitColls_ ){
      auto spColl = event.getValidHandle<StepPointMCCollection>(iColl);
      spMCColls.push_back( *spColl );
    }


    if(dumpSimParticleLeaves_) {
      auto ih = event.getValidHandle<SimParticleCollection>(input_);
      for(const auto& p : *ih) {
        if(is_leave(p.second)) {
          art::Ptr<SimParticle> pp(ih, p.first.asUint());
          process(pp, spMCColls);
        }
      }
    }
    else {
      auto ih = event.getValidHandle<SimParticlePtrCollection>(input_);
      for(const auto& p : *ih) {
        process(p, spMCColls);
      }
    }

  }

  //================================================================
  void StoppedParticlesDumper::process(const art::Ptr<SimParticle>& p, const VspMC& spMCColls) {
    const float tau = writeProperTime_ ? SimParticleGetTau::calculate(p,spMCColls,decayOffCodes_) : -1;
    data_ = StopInfo(p, spMCColls, tau);
    nt_->Fill();
  }

  //================================================================
  bool StoppedParticlesDumper::is_leave(const SimParticle& p) {
    return p.daughters().empty();
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::StoppedParticlesDumper);
