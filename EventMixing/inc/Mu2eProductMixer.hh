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
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/IO/ProductMix/MixHelper.h"

#include "Offline/MCDataProducts/inc/GenEventCount.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/MCTrajectoryCollection.hh"
#include "Offline/MCDataProducts/inc/CaloShowerStep.hh"
#include "Offline/MCDataProducts/inc/StrawGasStep.hh"
#include "Offline/MCDataProducts/inc/CrvStep.hh"
#include "Offline/MCDataProducts/inc/SurfaceStep.hh"
#include "Offline/MCDataProducts/inc/ExtMonFNALSimHit.hh"
#include "Offline/MCDataProducts/inc/CosmicLivetime.hh"
#include "Offline/DataProducts/inc/StageNormalization.hh"
#include "Offline/MCDataProducts/inc/SimTimeOffset.hh"
#include "Offline/MCDataProducts/inc/PhysicalVolumeInfoMultiCollection.hh"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"
#include "Offline/RecoDataProducts/inc/CrvDigi.hh"
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/RecoDataProducts/inc/CaloDigi.hh"
#include "Offline/MCDataProducts/inc/CaloShowerSim.hh"
#include "Offline/MCDataProducts/inc/CrvDigiMC.hh"
#include "Offline/DataProducts/inc/EventWindowMarker.hh"


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

    struct CosmicLivetimeMixerConfig {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag> moduleLabel{ Name("moduleLabel"), Comment("Input module label") };
      fhicl::Atom<std::string> srOutInstance{ Name("srOutInstance"), Comment("Output instance name for SubRun outputs"), "mixed" };
      fhicl::Atom<std::string> genCounterLabel{ Name("genCounterLabel"), Comment("Module label for the GenEventCounter"), "genCounter" };
    };

    // StageNormalization in SubRuns: carry the resampled pool's
    // origin-equivalent generated count through to the output, scaled by the
    // number of draws. Without this the chain breaks here -- the secondary's
    // GenEventCount does not reach the output, and the output's own records
    // draws, not generated events.
    struct StageNormMixerConfig {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      // EXACTLY ONE of moduleLabel or (poolGenCount + poolEventCount).
      // They are alternatives rather than a preference and a fallback because
      // art throws ProductNotFound when a declared mix op's product is absent
      // from the secondary -- it does not call the op with an empty input --
      // so declaring the StageNormalization op against a pool that has none
      // aborts the job on the first event.
      fhicl::OptionalAtom<art::InputTag> moduleLabel{ Name("moduleLabel"),
        Comment("SubRun StageNormalization of the RESAMPLED INPUT, as its "
                "StageNormalizationCounter wrote it. The pool MUST carry this "
                "product; use the genCounterLabel bootstrap if it does not.") };
      fhicl::Atom<std::string> srOutInstance{ Name("srOutInstance"),
        Comment("Output instance name for SubRun outputs"), "resampled" };
      // Bootstrap for a pool predating StageNormalization: state BOTH of the
      // pool's totals. Reading the generated count from the secondary instead
      // does not work -- the mix op sees the DRAWN SUBRUN's GenEventCount,
      // while the event count can only be stated for the pool as a whole, and
      // dividing one by the other mixes scopes and is wrong by the number of
      // subruns. Both totals are over ALL of fileNames; draws are uniform over
      // the pool, so the pool-level ratio is the correct per-draw expectation.
      fhicl::OptionalAtom<double> poolGenCount{ Name("poolGenCount"),
        Comment("Total GenEventCount over ALL of fileNames (the pool's generated "
                "events). Requires poolEventCount.") };
      fhicl::OptionalAtom<double> poolEventCount{ Name("poolEventCount"),
        Comment("Total events over ALL of fileNames. Requires poolGenCount.") };
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
      fhicl::Table<CollectionMixerConfig> surfaceStepMixer { fhicl::Name("surfaceStepMixer") };
      fhicl::Table<CollectionMixerConfig> extMonSimHitMixer { fhicl::Name("extMonSimHitMixer") };
      fhicl::Table<CollectionMixerConfig> eventIDMixer { fhicl::Name("eventIDMixer") };
      fhicl::Table<CollectionMixerConfig> strawDigiMixer { fhicl::Name("strawDigiMixer") };
      fhicl::Table<CollectionMixerConfig> strawDigiADCWaveformMixer { fhicl::Name("strawDigiADCWaveformMixer") };
      fhicl::Table<CollectionMixerConfig> strawDigiMCMixer { fhicl::Name("strawDigiMCMixer") };
      fhicl::Table<CollectionMixerConfig> caloDigiMixer { fhicl::Name("caloDigiMixer") };
      fhicl::Table<CollectionMixerConfig> caloShowerSimMixer { fhicl::Name("caloShowerSimMixer") };
      fhicl::Table<CollectionMixerConfig> crvDigiMixer { fhicl::Name("crvDigiMixer") };
      fhicl::Table<CollectionMixerConfig> crvDigiMCMixer { fhicl::Name("crvDigiMCMixer") };  //needs to be index-matched with crvDigiMixer
      fhicl::Table<CollectionMixerConfig> eventWindowMarkerMixer { fhicl::Name("eventWindowMarkerMixer") };
      fhicl::OptionalTable<CosmicLivetimeMixerConfig> cosmicLivetimeMixer { fhicl::Name("cosmicLivetimeMixer") };
      fhicl::OptionalTable<VolumeInfoMixerConfig> volumeInfoMixer { fhicl::Name("volumeInfoMixer") };
      fhicl::OptionalTable<StageNormMixerConfig> stageNormMixer { fhicl::Name("stageNormMixer") };
      fhicl::OptionalAtom<art::InputTag> simTimeOffset { fhicl::Name("simTimeOffset"), fhicl::Comment("Simulation time offset to apply (optional)") };
    };

    Mu2eProductMixer(const Config& conf, art::MixHelper& helper);

    void startEvent(art::Event const& e);
    void processEventIDs(const art::EventIDSequence& seq);
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

    bool mixSurfaceSteps(std::vector<SurfaceStepCollection const*> const& in,
                            SurfaceStepCollection& out,
                            art::PtrRemapper const& remap);

    bool mixExtMonSimHits(std::vector<ExtMonFNALSimHitCollection const*> const& in,
                          ExtMonFNALSimHitCollection& out,
                          art::PtrRemapper const& remap);

    bool mixStrawDigis(std::vector<StrawDigiCollection const*> const& in,
                       StrawDigiCollection& out,
                       art::PtrRemapper const& remap);

    bool mixStrawDigiADCWaveforms(std::vector<StrawDigiADCWaveformCollection const*> const& in,
                       StrawDigiADCWaveformCollection& out,
                       art::PtrRemapper const& remap);

    bool mixStrawDigiMCs(std::vector<StrawDigiMCCollection const*> const& in,
                       StrawDigiMCCollection& out,
                       art::PtrRemapper const& remap);

    bool mixCaloDigis(std::vector<CaloDigiCollection const*> const& in,
                       CaloDigiCollection& out,
                       art::PtrRemapper const& remap);

    bool mixCaloShowerSims(std::vector<CaloShowerSimCollection const*> const& in,
                       CaloShowerSimCollection& out,
                       art::PtrRemapper const& remap);

    bool mixCrvDigis(std::vector<CrvDigiCollection const*> const& in,
                       CrvDigiCollection& out,
                       art::PtrRemapper const& remap);

    bool mixCrvDigiMCs(std::vector<CrvDigiMCCollection const*> const& in,
                       CrvDigiMCCollection& out,
                       art::PtrRemapper const& remap);

    bool mixEventWindowMarkers(std::vector<EventWindowMarker const*> const& in,
                       EventWindowMarker& out,
                       art::PtrRemapper const& remap);


    bool mixEventIDs(std::vector<art::EventIDSequence const*> const &in,
                     art::EventIDSequence& out,
                     art::PtrRemapper const& remap);

    //----------------
    bool mixVolumeInfos(std::vector<PhysicalVolumeInfoMultiCollection const*> const& in,
                        PhysicalVolumeInfoMultiCollection& out,
                        art::PtrRemapper const&);

    bool mixGenEventCount(std::vector<GenEventCount const*> const& in,
                          GenEventCount& out,
                          art::PtrRemapper const& remap);

    bool mixCosmicLivetime(std::vector<mu2e::CosmicLivetime const*> const &in,
                           mu2e::CosmicLivetime& out,
                           art::PtrRemapper const& remap);

    bool mixStageNormalization(std::vector<mu2e::StageNormalization const*> const& in,
                               mu2e::StageNormalization& out,
                               art::PtrRemapper const& remap);

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

    typedef StrawGasStepCollection::size_type SGSOffset;
    std::vector<SGSOffset> sgsOffsets_;

    typedef CaloShowerStepCollection::size_type CSSOffset;
    std::vector<CSSOffset> cssOffsets_;
    typedef CrvStepCollection::size_type CSOffset;
    std::vector<CSOffset> csOffsets_;

    void updateSimParticle(SimParticle& particle, SPOffset offset, art::PtrRemapper const& remap);

    typedef std::map<cet::map_vector_key,PhysicalVolumeInfo> VolumeMap;
    typedef std::vector<VolumeMap> MultiStageMap;
    MultiStageMap subrunVolumes_;
    bool mixVolumes_;
    art::InputTag volumesInput_;
    std::string subrunVolInstanceName_;
    std::optional<std::string> evtVolInstanceName_;
    bool applyTimeOffset_;
    art::InputTag timeOffsetTag_;
    SimTimeOffset stoff_; // time offset for SimParticles and StepPointMCs

    void addInfo(VolumeMap* map, const PhysicalVolumeInfoSingleStage::value_type& entry);

    std::string subrunLivetimeInstanceName_;
    bool mixCosmicLivetimes_;
    unsigned int generatedEvents_ = 0;
    unsigned int resampledEvents_ = 0;
    unsigned int totalPrimaries_ = 0;
    float area_ = 0;
    float lowE_ = 0;
    float highE_ = 0;
    float fluxConstant_ = 0;
    float livetime_ = 0;
    // Merging cosmic livetime info from multiple subruns is not
    // implemented.  Make sure we only see a single subrun in a job.
    bool cosmicSubrunInitialized_ = false;
    art::SubRunID cosmicSubRun_;

    // StageNormalization mixing. perEventSum_ accumulates, once per drawn
    // event, the origin-generated events that draw stands for. Accumulating
    // rather than latching (as the cosmic path does) is what lets the draws
    // span many subruns of the pool, which they always do for a beam sample.
    std::string subrunStageNormInstanceName_;
    bool mixStageNorm_ = false;
    double perEventSum_ = 0.;
    uint64_t nDraws_ = 0;
    unsigned upstreamStages_ = 0;
    // Bootstrap for a pool with no StageNormalization of its own: both totals
    // are stated, and the draw count comes from resampledEvents_.
    bool bootstrapStageNorm_ = false;
    double poolGenCount_ = 0.;
    double poolEventCount_ = 0.;
    // Draws whose subrun carried no usable normalization: reported, never
    // quietly folded in as if they had one.
    uint64_t nDrawsUnnormalized_ = 0;

  };

}

#endif/*EventMixing_inc_Mu2eProductMixing_hh*/
