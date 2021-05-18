// A class that encapsulates event mixing code that is common for
// different Mu2e use cases.  It provides and registers callbacks to
// mix various data products, and is supposed to be used by MixFilter
// "detail" classes.  More documentation can be found in
//
//    art/Framework/Modules/MixFilter.h
//    art/Framework/IO/ProductMix/MixHelper.h
//    art/Framework/Core/PtrRemapper.h
//    art/Persistency/Common/CollectionUtilities.h
//
//
// Andrei Gaponenko, 2018

#ifndef EventMixing_inc_Mu2eProductMixing_hh
#define EventMixing_inc_Mu2eProductMixing_hh

#include <string>
#include <vector>
#include <optional>

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/OptionalTable.h"
#include "fhiclcpp/types/TupleAs.h"
#include "canvas/Utilities/InputTag.h"

#include "art/Framework/IO/ProductMix/MixHelper.h"

#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/MCTrajectoryCollection.hh"
#include "MCDataProducts/inc/CaloShowerStep.hh"
#include "MCDataProducts/inc/StrawGasStep.hh"
#include "MCDataProducts/inc/CrvStep.hh"
#include "MCDataProducts/inc/ExtMonFNALSimHitCollection.hh"
#include "MCDataProducts/inc/ProtonBunchIntensity.hh"
#include "MCDataProducts/inc/SimParticleTimeMap.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoMultiCollection.hh"



//================================================================
namespace mu2e {

  class Mu2eProductMixer {
  public:

    // Configuration for mixing one type of data products.
    struct CollectionMixerConfig {
      struct Entry {
        art::InputTag inTag;
        std::string outInstance;

        // Some outInstance inputs should not be treated literally,
        // this function is responsible for interpreting them.
        std::string resolvedInstanceName() const;

        Entry(const art::InputTag& i, const std::string& o): inTag(i), outInstance(o) {}
      };

      fhicl::Sequence<fhicl::TupleAs<Entry(art::InputTag,std::string)> >
      mixingMap { fhicl::Name("mixingMap"),
          fhicl::Comment("A sequence of InputTag to outputInstanceName"
                         " mappings for collections to be mixed."),
          std::vector<Entry>()
          };
    };

    // PhysicalVolumeInfoMultiCollection in SubRuns
    struct VolumeInfoMixerConfig {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag> srInput{ Name("srInput"), Comment("Input volume collection in SubRun") };
      fhicl::Atom<std::string> srOutInstance{ Name("srOutInstance"), Comment("Output instance name for SubRun outputs"), "" };

      fhicl::OptionalAtom<std::string> evtOutInstanceName { fhicl::Name("evtOutInstanceName"),
          Comment("If a value is provided, volume info  mixing will put a partial output into every event. ")
          };
    };

    // Configuration for the Mu2eProductMixing helper
    struct Config {
      fhicl::Table<CollectionMixerConfig> genParticleMixer { fhicl::Name("genParticleMixer") };
      fhicl::Table<CollectionMixerConfig> simParticleMixer { fhicl::Name("simParticleMixer") };
      fhicl::Table<CollectionMixerConfig> stepPointMCMixer { fhicl::Name("stepPointMCMixer") };
      fhicl::Table<CollectionMixerConfig> mcTrajectoryMixer { fhicl::Name("mcTrajectoryMixer") };
      fhicl::Table<CollectionMixerConfig> caloShowerStepMixer { fhicl::Name("caloShowerStepMixer") };
      fhicl::Table<CollectionMixerConfig> strawGasStepMixer { fhicl::Name("strawGasStepMixer") };
      fhicl::Table<CollectionMixerConfig> crvStepMixer { fhicl::Name("crvStepMixer") };
      fhicl::Table<CollectionMixerConfig> extMonSimHitMixer { fhicl::Name("extMonSimHitMixer") };
      fhicl::Table<CollectionMixerConfig> protonBunchIntensityMixer { fhicl::Name("protonBunchIntensityMixer") };
      fhicl::Table<CollectionMixerConfig> protonTimeMapMixer { fhicl::Name("protonTimeMapMixer") };
      fhicl::Table<CollectionMixerConfig> eventIDMixer { fhicl::Name("eventIDMixer") };

      fhicl::OptionalTable<VolumeInfoMixerConfig> volumeInfoMixer { fhicl::Name("volumeInfoMixer") };
    };

    Mu2eProductMixer(const Config& conf, art::MixHelper& helper);

    void beginSubRun(const art::SubRun& sr);
    void endSubRun(art::SubRun& sr);

  private:

    bool mixGenParticles(std::vector<GenParticleCollection const*> const& in,
                         GenParticleCollection& out,
                         art::PtrRemapper const& remap);

    bool mixSimParticles(std::vector<SimParticleCollection const*> const& in,
                         SimParticleCollection& out,
                         art::PtrRemapper const& remap);

    bool mixStepPointMCs(std::vector<StepPointMCCollection const*> const& in,
                         StepPointMCCollection& out,
                         art::PtrRemapper const& remap);

    bool mixMCTrajectories(std::vector<MCTrajectoryCollection const*> const& in,
                           MCTrajectoryCollection& out,
                           art::PtrRemapper const& remap);

    bool mixCaloShowerSteps(std::vector<CaloShowerStepCollection const*> const& in,
                            CaloShowerStepCollection& out,
                            art::PtrRemapper const& remap);

    bool mixStrawGasSteps(std::vector<StrawGasStepCollection const*> const& in,
                            StrawGasStepCollection& out,
                            art::PtrRemapper const& remap);

    bool mixCrvSteps(std::vector<CrvStepCollection const*> const& in,
                            CrvStepCollection& out,
                            art::PtrRemapper const& remap);

    bool mixExtMonSimHits(std::vector<ExtMonFNALSimHitCollection const*> const& in,
                          ExtMonFNALSimHitCollection& out,
                          art::PtrRemapper const& remap);

    bool mixProtonBunchIntensity(std::vector<mu2e::ProtonBunchIntensity const*> const &in,
                                 mu2e::ProtonBunchIntensity& out,
                                 art::PtrRemapper const& remap);

    bool mixProtonTimeMap(std::vector<mu2e::SimParticleTimeMap const*> const &in,
                                mu2e::SimParticleTimeMap& out,
                                art::PtrRemapper const& remap);

    bool mixEventIDs(std::vector<art::EventIDSequence const*> const &in,
                     art::EventIDSequence& out,
                     art::PtrRemapper const& remap);


    //----------------
    bool mixVolumeInfos(std::vector<PhysicalVolumeInfoMultiCollection const*> const& in,
                        PhysicalVolumeInfoMultiCollection& out,
                        art::PtrRemapper const&);

    //----------------
    // If elements of a collection can be pointed to by other
    // collections, the offset array for the pointed-to collection
    // need to be a data member here, as the offsets will used by more
    // than one mixOp.
    typedef SimParticleCollection::size_type SPOffset;
    typedef std::vector<SPOffset> SPOffsets;
    SPOffsets simOffsets_;

    typedef GenParticleCollection::size_type GenOffset;
    std::vector<GenOffset> genOffsets_;

    void updateSimParticle(SimParticle& particle, SPOffset offset, art::PtrRemapper const& remap);

    typedef std::map<cet::map_vector_key,PhysicalVolumeInfo> VolumeMap;
    typedef std::vector<VolumeMap> MultiStageMap;
    MultiStageMap subrunVolumes_;
    bool mixVolumes_;
    art::InputTag volumesInput_;
    std::string subrunVolInstanceName_;
    std::optional<std::string> evtVolInstanceName_;

    void addInfo(VolumeMap* map, const PhysicalVolumeInfoSingleStage::value_type& entry);

  };

}

#endif/*EventMixing_inc_Mu2eProductMixing_hh*/
