// Andrei Gaponenko, 2015, 2021

#ifndef Mu2eG4_Mu2eG4Inputs_hh
#define Mu2eG4_Mu2eG4Inputs_hh

#include <string>
#include <vector>

#include "canvas/Utilities/InputTag.h"
#include "Mu2eG4/inc/Mu2eG4Config.hh"
#include "Mu2eG4/inc/Mu2eG4PrimaryType.hh"

#include "art/Framework/Principal/Handle.h"
#include "MCDataProducts/inc/SimParticleCollection.hh"

namespace art { class Event; }

namespace mu2e {

  class Mu2eG4Inputs {
  public:
    explicit Mu2eG4Inputs(const Mu2eG4Config::Inputs_& conf);

    bool multiStage() const { return multiStage_; }

    Mu2eG4PrimaryType primaryType() const { return primaryType_; }

    const art::InputTag& primaryTag() const { return primaryTag_; }

    const art::InputTag& inputMCTrajectories() const { return inputMCTrajectories_; }

    const art::InputTag& inputPhysVolumeMultiInfo() const { return inputPhysVolumeMultiInfo_; }

    // The handle is not valid if there are no input sim particles, either for
    // initial stage jobs, or for non-filtered multistage input when the
    // input event has no primaries for the current stage.
    art::Handle<SimParticleCollection> inputSimParticles(const art::Event& evt) const;

  private:
    Mu2eG4PrimaryType primaryType_;
    art::InputTag primaryTag_;
    art::InputTag inputMCTrajectories_;
    art::InputTag inputPhysVolumeMultiInfo_;

    // derived
    bool multiStage_;
  };

} // end namespace mu2e

#endif /* Mu2eG4_Mu2eG4Inputs_hh */
