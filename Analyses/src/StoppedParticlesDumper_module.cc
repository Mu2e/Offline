// Write stopped particles into an ntuple.
//
// Andrei Gaponenko, 2013

#include <string>
#include <algorithm>

#include "cetlib/exception.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Utilities/InputTag.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/SimParticlePtrCollection.hh"
#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"

#include "ConditionsService/inc/GlobalConstantsHandle.hh"
#include "ConditionsService/inc/PhysicsParams.hh"

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
    art::InputTag input_;
    bool writeProperTime_;
    VS hitColls_;

    std::vector<int> decayOffCodes_;

    GlobalConstantsHandle<PhysicsParams> gc_;

    TTree *nt_;
    StopInfo data_;

    float getTau(const art::Ptr<SimParticle>& p, const VspMC& hitColls);
    bool decayTurnedOff(PDGCode::type pdgId);
  };

  //================================================================
  StoppedParticlesDumper::StoppedParticlesDumper(const fhicl::ParameterSet& pset) :
    art::EDAnalyzer(pset),
    input_(pset.get<std::string>("inputCollection")),
    writeProperTime_(pset.get<bool>("writeProperTime")),
    hitColls_( writeProperTime_ ? pset.get<VS>("hitCollections") : VS() ),
    nt_()
  {
    if(writeProperTime_) {
      decayOffCodes_ = pset.get<std::vector<int> >("decayOffPDGCodes");
      // must sort to use binary_search
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

    auto ih = event.getValidHandle<SimParticlePtrCollection>(input_);
    for(const auto& p : *ih) {
      const float tau = writeProperTime_ ? getTau(p,spMCColls) : -1;
      data_ = StopInfo( p, spMCColls, tau);
      nt_->Fill();
    }

  }

  //================================================================
  float StoppedParticlesDumper::getTau( const art::Ptr<SimParticle>& p,
                                        const VspMC& hitColls ){

    double tau = p->endProperTime() / gc_->getParticleLifetime(p->pdgId());;

    // The mu2ePrimary code means that G4 track was created by our PrimaryGeneratorAction.
    // If the particle is "mu2ePrimary" but still has a parent, it is a continuation
    // of a particle from the previous simulation stage, and their proper times should
    // be combined.

    art::Ptr<SimParticle> part (p);
    while(part->parent().isNonnull()) {

      if((part->creationCode() == ProcessCode::mu2ePrimary)) {

        // The current particle is a continuation from the previous stage,
        // not a physically different particle.  We need to find its record
        // in the StepPointMC collections to get the correct proper time.
        part = part->parent();

        // Find matches in hit collections
        unsigned counter (0);
        const StepPointMC* spMC(nullptr);
        for ( const auto& hitColl : hitColls ) {
          std::for_each( hitColl.begin(), hitColl.end(),
                         [&](const StepPointMC& sp){
                           if ( sp.simParticle().key() == part.key() ) {
                             spMC = &sp;
                             counter++;
                           }
                         } );

        }

        if      ( counter == 0 ) throw cet::exception("StepPointMC") << " Non-existent StepPointMC-SimParticle assignment! " ;
        else if ( counter  > 1 ) throw cet::exception("StepPointMC") << " Ambiguous StepPointMC-SimParticle assignment! " ;
        else  tau += spMC->properTime() / gc_->getParticleLifetime(part->pdgId());;
      }
      else {
        // The current particle was produced by a G4 physics process.
        // See if proper time of its ancestor should be included.
        part = part->parent();
        if(decayTurnedOff(part->pdgId())) {
          tau += part->endProperTime() / gc_->getParticleLifetime(part->pdgId());
        }
      }
    } // loop up to the primary

    return tau;
  }

  //================================================================
  bool StoppedParticlesDumper::decayTurnedOff(PDGCode::type pdgId) {
    return std::binary_search(decayOffCodes_.begin(), decayOffCodes_.end(), int(pdgId));
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::StoppedParticlesDumper);
