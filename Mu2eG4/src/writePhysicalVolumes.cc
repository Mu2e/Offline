// Andrei Gaponenko, 2021

#include "Offline/Mu2eG4/inc/writePhysicalVolumes.hh"

#include <memory>

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"

namespace {

  // Helper functions to work around the problem that, starting in art v3.11,
  // the signature of put is different for events and subruns.
  // A better solution would be to use template logic to make the decision.
  using Collection_t = mu2e::PhysicalVolumeInfoMultiCollection;

  void storeIt  (art::SubRun& sr, std::unique_ptr<Collection_t>& up, const std::string& outInstanceName ){
    sr.put( std::move(up), outInstanceName, art::fullSubRun());
  }

  void storeIt  (art::Event& e, std::unique_ptr<Collection_t>& up, const std::string& outInstanceName ){
    e.put( std::move(up), outInstanceName);
  }

}

namespace mu2e {

  template <class PRINCIPAL>
  unsigned writePhysicalVolumes(PRINCIPAL& store,
                                const std::optional<art::InputTag>& input,
                                const PhysicalVolumeInfoSingleStage& vi,
                                const std::string& outInstanceName)
  {
    using Collection_t = PhysicalVolumeInfoMultiCollection;
    auto mvi = std::make_unique<Collection_t>();

    if(input) {
      // Copy over data from the previous simulation stages
      auto const& ih = store.template getValidHandle<Collection_t>(*input);
      mvi->reserve(1 + ih->size());
      mvi->insert(mvi->begin(), ih->cbegin(), ih->cend());
    }

    // By definition simStage=0 if we start with GenParticles and not doing multiStage.
    unsigned const simStage = mvi->size();

    // Append info for the current stage
    mvi->emplace_back(vi);

    storeIt ( store, mvi, outInstanceName);

    return simStage;
  }

  //----------------------------------------------------------------
  // Use explicit instantiation instead of puting function definition into a header

  template unsigned
  writePhysicalVolumes<art::SubRun>(art::SubRun&,
                                    const std::optional<art::InputTag>& input,
                                    const PhysicalVolumeInfoSingleStage& vi,
                                    const std::string& outInstanceName);

  template unsigned
  writePhysicalVolumes<art::Event>(art::Event&,
                                   const std::optional<art::InputTag>& input,
                                   const PhysicalVolumeInfoSingleStage& vi,
                                   const std::string& outInstanceName);

  //----------------------------------------------------------------
}
