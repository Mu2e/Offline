// Andrei Gaponenko, 2015, 2021

#ifndef Mu2eG4_Mu2eG4Inputs_hh
#define Mu2eG4_Mu2eG4Inputs_hh

#include <string>
#include <vector>
#include <optional>

#include "cetlib/maybe_ref.h"
#include "canvas/Utilities/InputTag.h"
#include "Offline/Mu2eG4/inc/Mu2eG4Config.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4PrimaryType.hh"

#include "art/Framework/Principal/Handle.h"
#include "Offline/MCDataProducts/inc/SimParticle.hh"

namespace art { class Event; }

namespace mu2e {

  struct EventLevelVolInfos {
    art::InputTag input;
    std::string outInstance;
  };

  class Mu2eG4Inputs {
  public:

    // The two elements of the return value of the function inputSimParticles.
    // The maybe_ref is not valid if there are no input sim particles, either for
    // initial stage jobs, or for non-filtered multistage input when the
    // input event has no primaries for the current stage.
    struct InputSimsInfo{
      cet::maybe_ref<cet::map_vector<mu2e::SimParticle> const > sims;
      art::ProductID                                            id;

      InputSimsInfo ():
        sims(), id(){
      }

      InputSimsInfo ( cet::map_vector<mu2e::SimParticle> const&  asims
                      , art::ProductID const& aid):
        sims(asims), id(aid){
      }

      void reseat( cet::map_vector<mu2e::SimParticle> const&  asims
                   , art::ProductID const& aid){
        sims.reseat(asims);
        id   = aid;
      }

      bool isValid() const { return sims.isValid(); };
    };

    explicit Mu2eG4Inputs(const Mu2eG4Config::Inputs_& conf);

    bool multiStage() const { return multiStage_; }

    Mu2eG4PrimaryType primaryType() const { return primaryType_; }

    const art::InputTag& primaryTag() const { return primaryTag_; }

    const art::InputTag& inputMCTrajectories() const { return inputMCTrajectories_; }

    std::optional<unsigned> simStageOverride() const { return simStageOverride_; }

    const art::InputTag& inputPhysVolumeMultiInfo() const { return inputPhysVolumeMultiInfo_; }

    const std::optional<EventLevelVolInfos>& updateEventLevelVolumeInfos() const { return elvi_; }

    InputSimsInfo inputSimParticles(const art::Event& evt) const;

  private:
    Mu2eG4PrimaryType primaryType_;
    art::InputTag primaryTag_;
    art::InputTag inputMCTrajectories_;
    std::optional<unsigned> simStageOverride_;
    art::InputTag inputPhysVolumeMultiInfo_;
    std::optional<EventLevelVolInfos> elvi_;
    // derived
    bool multiStage_;
  };

} // end namespace mu2e

#endif /* Mu2eG4_Mu2eG4Inputs_hh */
