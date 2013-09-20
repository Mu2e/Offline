// Write stopped particles into an ntuple.
//
// Andrei Gaponenko, 2013

#include <string>

#include "cetlib/exception.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Utilities/InputTag.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/SimParticlePtrCollection.hh"

#include "TTree.h"

namespace mu2e {

  namespace {
    // This should be minimal info, we'll want to load this
    // in memory in consumer jobs.  This is NOT an analysis ntuple!
    struct StopInfo {
      float x;
      float y;
      float z;
      float t;
      float tau; // proper time, for stopped pion weights

      StopInfo() : x(), y(), z(), t(), tau() {}

      StopInfo(const SimParticle& p)
        : x(p.endPosition().x())
        , y(p.endPosition().y())
        , z(p.endPosition().z())
        , t(p.endGlobalTime())
        , tau(p.endProperTime())
      {
        if(!p.endDefined()) {
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
    art::InputTag input_;
    bool writeProperTime_;

    TTree *nt_;
    StopInfo data_;
  };

  //================================================================
  StoppedParticlesDumper::StoppedParticlesDumper(const fhicl::ParameterSet& pset)
    : input_(pset.get<std::string>("inputCollection"))
    , writeProperTime_(pset.get<bool>("writeProperTime"))
    , nt_()
  {}

  //================================================================
  void StoppedParticlesDumper::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;
    std::string branchDesc("x/F:y/F:z/F:time/F");
    if(writeProperTime_) {
      branchDesc += ":tau/F";
    }
    nt_ = tfs->make<TTree>( "stops", "Stopped particles ntuple");
    nt_->Branch("stops", &data_, branchDesc.c_str());
  }

  //================================================================
  void StoppedParticlesDumper::analyze(const art::Event& event) {
    auto ih = event.getValidHandle<SimParticlePtrCollection>(input_);
    for(const auto& p : *ih) {
      data_ = StopInfo(*p);
      nt_->Fill();
    }
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::StoppedParticlesDumper);
