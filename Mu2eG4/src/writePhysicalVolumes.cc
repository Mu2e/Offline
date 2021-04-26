// Andrei Gaponenko, 2021

#include "Mu2eG4/inc/writePhysicalVolumes.hh"

#include <memory>

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"

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

    store.put(std::move(mvi), outInstanceName);

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
