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
#include "art_root_io/TFileService.h"

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

    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::Atom<bool> dumpSimParticleLeaves {
        Name("dumpSimParticleLeaves"),
          Comment("The mode selector.  Set it to true to dump info on all the leaves\n"
                  "of a SimParticleCollection.   If it is set to false, the input should\n"
                  "be of the SimParticlePtrCollection type, and its full content will be included.\n"
                  ),
          false
          };

      fhicl::Atom<art::InputTag> inputCollection {
        Name("inputCollection"),
          Comment("InputTag of the collection to use.")
          };

      fhicl::Atom<bool> writeProperTime {
        Name("writeProperTime"),
          Comment("Compute and write out normalized proper time tau of particles.\n"
                  "This quantity is defined in such a way that exp(-tau) is the survival\n"
                  "probability of the particle.  It is used in cases when the decay\n"
                  "process is disabled during the simulation."
                  ),
          false
          };

      fhicl::Sequence<int> decayOffPDGCodes {
        Name("decayOffPDGCodes"),
          Comment("A list of PDG IDs of particles that had their decay process turned off during\n"
                  "the simulation. This must be specified if and only if writeProperTime is requested."),
          [this](){ return writeProperTime(); }
      };

      fhicl::Sequence<art::InputTag> hitCollections {
        Name("hitCollections"),
          Comment("A list of StepPointMCCollection-s via which different stages of simulation\n"
                  "are connected.  Must be specified if and only if writeProperTime is requested."),
          [this](){ return writeProperTime(); }
      };

    };

    using Parameters = art::EDAnalyzer::Table<Config>;
    explicit StoppedParticlesDumper(const Parameters& conf);

    void beginJob() override;
    void analyze(const art::Event& evt) override;
  private:
    bool dumpSimParticleLeaves_;
    art::InputTag input_;
    bool writeProperTime_;
    std::vector<art::InputTag> hitColls_;

    std::vector<int> decayOffCodes_;

    TTree *nt_;
    StopInfo data_;

    bool is_leave(const SimParticle& p);
    void process(const art::Ptr<SimParticle>& p, const VspMC& spMCcolls);
  };

  //================================================================
  StoppedParticlesDumper::StoppedParticlesDumper(const Parameters& conf) :
    art::EDAnalyzer(conf),
    dumpSimParticleLeaves_(conf().dumpSimParticleLeaves()),
    input_(conf().inputCollection()),
    writeProperTime_(conf().writeProperTime()),
    nt_()
  {
    if(writeProperTime_) {
      hitColls_ =  conf().hitCollections();
      decayOffCodes_ = conf().decayOffPDGCodes();

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
