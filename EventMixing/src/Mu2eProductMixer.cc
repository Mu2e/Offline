// Andrei Gaponenko, 2018

#include "Offline/EventMixing/inc/Mu2eProductMixer.hh"

#include <utility>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <string>

#include "cetlib_except/exception.h"

#include "art/Persistency/Common/CollectionUtilities.h"

//================================================================
namespace mu2e {

  //----------------------------------------------------------------
  std::string Mu2eProductMixer::CollectionMixerConfig::Entry::resolvedInstanceName() const {
    return (outInstance == ":") ? inTag.instance() :  outInstance;
  }

  //----------------------------------------------------------------
  namespace {

    // newEntryIndex is an address in the output (flattened)
    // collection and the offsets are those that were recorded by the
    // flattenCollections() call
    template<typename IndexType, typename OFFSETS>
    typename OFFSETS::size_type
    getInputEventIndex(IndexType newEntryIndex, const OFFSETS& offsets) {
      auto ub = std::upper_bound(offsets.begin(), offsets.end(), newEntryIndex);
      if(ub == offsets.begin()) {
        throw cet::exception("RANGE")<<"getInputEventIndex(): newEntryIndex="
                                     <<newEntryIndex
                                     <<" is below the first offset="
                                     <<*offsets.begin()
                                     <<std::endl;
      }
      return std::distance(offsets.begin(), --ub);
    }
  }

  //----------------------------------------------------------------
  Mu2eProductMixer::Mu2eProductMixer(const Config& conf, art::MixHelper& helper)
    : mixVolumes_(false)
      , applyTimeOffset_(conf.simTimeOffset.hasValue())
      , stoff_(0.0)
      , mixCosmicLivetimes_(false)
  {
    if(applyTimeOffset_){
      timeOffsetTag_ = conf.simTimeOffset().value();
    }

    for(const auto& e: conf.genParticleMixer().mixingMap()) {
      helper.declareMixOp
              (e.inTag, e.resolvedInstanceName(), &Mu2eProductMixer::mixGenParticles, *this);
    }

    for(const auto& e: conf.simParticleMixer().mixingMap()) {
      helper.declareMixOp
        (e.inTag, e.resolvedInstanceName(), &Mu2eProductMixer::mixSimParticles, *this);
    }

    for(const auto& e: conf.stepPointMCMixer().mixingMap()) {
      helper.declareMixOp
        (e.inTag, e.resolvedInstanceName(), &Mu2eProductMixer::mixStepPointMCs, *this);
    }

    for(const auto& e: conf.mcTrajectoryMixer().mixingMap()) {
      helper.declareMixOp
        (e.inTag, e.resolvedInstanceName(), &Mu2eProductMixer::mixMCTrajectories, *this);
    }

    for(const auto& e: conf.caloShowerStepMixer().mixingMap()) {
      helper.declareMixOp
        (e.inTag, e.resolvedInstanceName(), &Mu2eProductMixer::mixCaloShowerSteps, *this);
    }

    for(const auto& e: conf.strawGasStepMixer().mixingMap()) {
      helper.declareMixOp
        (e.inTag, e.resolvedInstanceName(), &Mu2eProductMixer::mixStrawGasSteps, *this);
    }

    for(const auto& e: conf.crvStepMixer().mixingMap()) {
      helper.declareMixOp
        (e.inTag, e.resolvedInstanceName(), &Mu2eProductMixer::mixCrvSteps, *this);
    }

    for(const auto& e: conf.extMonSimHitMixer().mixingMap()) {
      helper.declareMixOp
        (e.inTag, e.resolvedInstanceName(), &Mu2eProductMixer::mixExtMonSimHits, *this);
    }

    for(const auto& e: conf.eventIDMixer().mixingMap()) {
      helper.declareMixOp
        (e.inTag, e.resolvedInstanceName(), &Mu2eProductMixer::mixEventIDs, *this);
    }

    for(const auto& e: conf.strawDigiMixer().mixingMap()) {
      helper.declareMixOp
        (e.inTag, e.resolvedInstanceName(), &Mu2eProductMixer::mixStrawDigis, *this);
    }

    for(const auto& e: conf.strawDigiADCWaveformMixer().mixingMap()) {
      helper.declareMixOp
        (e.inTag, e.resolvedInstanceName(), &Mu2eProductMixer::mixStrawDigiADCWaveforms, *this);
    }

    for(const auto& e: conf.strawDigiMCMixer().mixingMap()) {
      helper.declareMixOp
        (e.inTag, e.resolvedInstanceName(), &Mu2eProductMixer::mixStrawDigiMCs, *this);
    }

    for(const auto& e: conf.eventWindowMarkerMixer().mixingMap()) {
      helper.declareMixOp
        (e.inTag, e.resolvedInstanceName(), &Mu2eProductMixer::mixEventWindowMarkers, *this);
    }

    //----------------------------------------------------------------
    // VolumeInfo handling

    VolumeInfoMixerConfig vmc;
    if(conf.volumeInfoMixer(vmc)) {

      mixVolumes_ = true;
      volumesInput_ = vmc.srInput();
      subrunVolInstanceName_ = vmc.srOutInstance();


      helper.produces<PhysicalVolumeInfoMultiCollection, art::InSubRun>(subrunVolInstanceName_);

      std::string tmp;
      if(vmc.evtOutInstanceName(tmp)) {
        evtVolInstanceName_.emplace(tmp);
      }

      const bool putVolsIntoEvent{evtVolInstanceName_};
      std::string evtOutInstance = putVolsIntoEvent ? *evtVolInstanceName_ : "unused";
      helper.declareMixOp<art::InSubRun>
        (volumesInput_, evtOutInstance, &Mu2eProductMixer::mixVolumeInfos, *this, putVolsIntoEvent);
    }

    //----------------------------------------------------------------
    // CosmicLivetime handling
    CosmicLivetimeMixerConfig clmc;
    if (conf.cosmicLivetimeMixer(clmc)) {
      mixCosmicLivetimes_ = true;
      subrunLivetimeInstanceName_ = clmc.srOutInstance();
      helper.produces<CosmicLivetime, art::InSubRun>(subrunLivetimeInstanceName_);
      helper.declareMixOp<art::InSubRun>
        (clmc.genCounterLabel(), "", &Mu2eProductMixer::mixGenEventCount, *this, false);
      helper.declareMixOp<art::InSubRun>
        (clmc.moduleLabel(), "", &Mu2eProductMixer::mixCosmicLivetime, *this);
    }

  }

  //================================================================
  void Mu2eProductMixer::startEvent(art::Event const& e) {
    if(applyTimeOffset_){
    // find the time offset in the event, and copy it locally
      const auto& stoH = e.getValidHandle<SimTimeOffset>(timeOffsetTag_);
      stoff_ = *stoH;
    }
    resampledEvents_++;
  }

  //----------------------------------------------------------------
  void Mu2eProductMixer::processEventIDs(const art::EventIDSequence& seq)  {
    if(mixCosmicLivetimes_) {

      if(seq.size() != 1) {
        throw cet::exception("BADINPUT")<<"Mu2eProductMixer: can't mix CosmicLiveTime" << std::endl;
      }

      if(!cosmicSubrunInitialized_) {
        cosmicSubrunInitialized_ = true;
        cosmicSubRun_ = seq.at(0).subRunID();
      }
      else {
        if(seq.at(0).subRunID() != cosmicSubRun_) {
          throw cet::exception("BADINPUT")<<"Mu2eProductMixer: input from multiple subruns is not supported for CosmicLivetime. "
                                          <<"Got SubRunIDs"<<cosmicSubRun_<<" and "<<seq.at(0).subRunID()
                                          <<std::endl;

        }
      }

    }
  }

  //----------------------------------------------------------------
  void Mu2eProductMixer::beginSubRun(const art::SubRun& s) {
    subrunVolumes_.clear();
    resampledEvents_ = 0;

  }

  //----------------------------------------------------------------
  void Mu2eProductMixer::endSubRun(art::SubRun& sr) {
    if(mixVolumes_) {
      auto col =  std::make_unique<PhysicalVolumeInfoMultiCollection>();

      col->resize(subrunVolumes_.size());
      for(unsigned stage=0; stage<subrunVolumes_.size(); ++stage) {
        (*col)[stage].insert(subrunVolumes_[stage].begin(), subrunVolumes_[stage].end());
      }

      sr.put(std::move(col), subrunVolInstanceName_, art::fullSubRun());
    }
    if (mixCosmicLivetimes_) {
      if(generatedEvents_ == 0)throw cet::exception("BADINPUT")<<"Mu2eProductMixer: generated event count =0; was the mixin file opened correctly?" << std::endl;
      float scaling = resampledEvents_ / generatedEvents_;
      auto livetime = std::make_unique<CosmicLivetime>(totalPrimaries_ * scaling,
                                                       area_, lowE_, highE_, fluxConstant_, livetime_ * scaling);
      sr.put(std::move(livetime), subrunLivetimeInstanceName_, art::fullSubRun());
    }
  }

  //----------------------------------------------------------------
  bool Mu2eProductMixer::mixGenParticles(std::vector<GenParticleCollection const*> const& in,
                                         GenParticleCollection& out,
                                         art::PtrRemapper const& remap)
  {
    art::flattenCollections(in, out, genOffsets_);
    if(applyTimeOffset_){
      for(auto& particle : out){
        particle.time() += stoff_.timeOffset_;
        // proper times are WRT the particles own internal clock, can't shift them
      }
    }


    return true;
  }

  //----------------------------------------------------------------
  bool Mu2eProductMixer::mixSimParticles(std::vector<SimParticleCollection const*> const& in,
                                         SimParticleCollection& out,
                                         art::PtrRemapper const& remap)
  {
    art::flattenCollections(in, out, simOffsets_ );

    // Update the Ptrs inside each SimParticle
    for(auto& entry: out) {
      auto key = entry.first;
      auto& particle = entry.second;
      auto ie = getInputEventIndex(key.asUint(), simOffsets_);
      updateSimParticle(particle, ie, remap);
    }
    return true;
  }

  //----------------
  // Update one SimParticle to deal with the flattening of the SimParticleCollections.
  void Mu2eProductMixer::updateSimParticle(mu2e::SimParticle& sim,
                                           SPOffsets::size_type inputEventIndex,
                                           art::PtrRemapper const& remap
                                           )
  {
    auto simOffset = simOffsets_[inputEventIndex];

    // Id of the SimParticle - must match the map_vector key.
    sim.id() = SimParticle::key_type( sim.id().asInt() + simOffset );

    // Ptr to the parent SimParticle.
    if ( sim.parent().isNonnull() ){
      sim.parent() = remap(sim.parent(), simOffset);
    }

    // Ptrs to all of the daughters.
    for(auto& d: sim.daughters()) {
      d = remap(d, simOffset);
    }

    // If we mix GenParticles, update that Ptr, too.
    if(!genOffsets_.empty()) {
      sim.genParticle() = remap( sim.genParticle(), genOffsets_[inputEventIndex]);
    }

    if(applyTimeOffset_){
      sim.startGlobalTime() += stoff_.timeOffset_;
      sim.endGlobalTime() += stoff_.timeOffset_;
      // proper times are WRT the particles own internal clock, can't shift them
    }

  }

  //----------------------------------------------------------------
  bool Mu2eProductMixer::mixStepPointMCs(std::vector<StepPointMCCollection const*> const& in,
                                         StepPointMCCollection& out,
                                         art::PtrRemapper const& remap)
  {
    std::vector<StepPointMCCollection::size_type> stepOffsets;
    art::flattenCollections(in, out, stepOffsets);

    for(StepPointMCCollection::size_type i=0; i<out.size(); ++i) {
      auto ie = getInputEventIndex(i, stepOffsets);
      auto& step = out[i];
      step.simParticle() = remap(step.simParticle(), simOffsets_[ie]);
      if(applyTimeOffset_){
        step.time() += stoff_.timeOffset_;
      }
    }
    return true;
  }

  //----------------------------------------------------------------
  bool Mu2eProductMixer::mixMCTrajectories(std::vector<MCTrajectoryCollection const*> const& in,
                                           MCTrajectoryCollection& out,
                                           art::PtrRemapper const& remap)
  {
    // flattenCollections() does not seem to preserve enough info to remap ptrs in the output map.
    // Follow the pattern, including the nullptr checks, but add custom remapping code
    std::pair<MCTrajectoryCollection::iterator,bool> res;
    for(std::vector<MCTrajectoryCollection const*>::size_type ieIndex = 0; ieIndex < in.size(); ++ieIndex) {
      if (in[ieIndex] != nullptr) {
        for(const auto & orig : *in[ieIndex]) {
          if(!applyTimeOffset_) {
            res = out.insert(std::make_pair(remap(orig.first, simOffsets_[ieIndex]),
                  MCTrajectory(remap(orig.second.sim(), simOffsets_[ieIndex]), orig.second.points())));
          } else {
            // make a deep copy of the points with shifted time
            std::vector<MCTrajectoryPoint> newpoints;
            newpoints.reserve(orig.second.points().size());
            for(auto const& mcpt : orig.second.points())
              newpoints.emplace_back(mcpt.pos(),mcpt.t()+stoff_.timeOffset_,mcpt.kineticEnergy());
            res = out.insert(std::make_pair(remap(orig.first, simOffsets_[ieIndex]),
                  MCTrajectory(remap(orig.second.sim(), simOffsets_[ieIndex]), newpoints)));
          }
          if(!res.second) {
            throw cet::exception("BUG")<<"mixMCTrajectories(): failed to insert an entry, ieIndex="<<ieIndex
              <<", orig ptr = "<<orig.first
              <<std::endl;
          }
        }
      }
    }

    return true;
  }

  //----------------------------------------------------------------
  bool Mu2eProductMixer::mixCaloShowerSteps(std::vector<CaloShowerStepCollection const*> const& in,
                                            CaloShowerStepCollection& out,
                                            art::PtrRemapper const& remap)
  {
    std::vector<CaloShowerStepCollection::size_type> stepOffsets;
    art::flattenCollections(in, out, stepOffsets);

    for(CaloShowerStepCollection::size_type i=0; i<out.size(); ++i) {
      auto ie = getInputEventIndex(i, stepOffsets);
      auto& step = out[i];
      step.setSimParticle( remap(step.simParticle(), simOffsets_[ie]) );
      if(applyTimeOffset_){
        step.time() += stoff_.timeOffset_;
      }
    }

    return true;
  }

  //----------------------------------------------------------------
  bool Mu2eProductMixer::mixStrawGasSteps(std::vector<StrawGasStepCollection const*> const& in,
                                          StrawGasStepCollection& out,
                                          art::PtrRemapper const& remap)
  {
    std::vector<StrawGasStepCollection::size_type> stepOffsets;
    art::flattenCollections(in, out, stepOffsets);

    for(StrawGasStepCollection::size_type i=0; i<out.size(); ++i) {
      auto ie = getInputEventIndex(i, stepOffsets);
      auto& step = out[i];
      step.simParticle() = remap(step.simParticle(), simOffsets_[ie]);
      if(applyTimeOffset_){
        step.time() += stoff_.timeOffset_;
      }
    }

    return true;
  }

  bool Mu2eProductMixer::mixCrvSteps(std::vector<CrvStepCollection const*> const& in,
                                          CrvStepCollection& out,
                                          art::PtrRemapper const& remap)
  {
    std::vector<CrvStepCollection::size_type> stepOffsets;
    art::flattenCollections(in, out, stepOffsets);

    for(CrvStepCollection::size_type i=0; i<out.size(); ++i) {
      auto ie = getInputEventIndex(i, stepOffsets);
      auto& step = out[i];
      step.simParticle() = remap(step.simParticle(), simOffsets_[ie]);
      if(applyTimeOffset_){
        step.startTime() += stoff_.timeOffset_;
        step.endTime() += stoff_.timeOffset_;
      }
    }

    return true;
  }

  //----------------------------------------------------------------
  bool Mu2eProductMixer::mixExtMonSimHits(std::vector<ExtMonFNALSimHitCollection const*> const& in,
                                          ExtMonFNALSimHitCollection& out,
                                          art::PtrRemapper const& remap)
  {
    std::vector<ExtMonFNALSimHitCollection::size_type> stepOffsets;
    art::flattenCollections(in, out, stepOffsets);

    for(ExtMonFNALSimHitCollection::size_type i=0; i<out.size(); ++i) {
      auto ie = getInputEventIndex(i, stepOffsets);
      auto& step = out[i];
      step.setSimParticle( remap(step.simParticle(), simOffsets_[ie]) );
    }

    return true;
  }

  bool Mu2eProductMixer::mixStrawDigis(std::vector<StrawDigiCollection const*> const& in,
                     StrawDigiCollection& out,
                     art::PtrRemapper const& remap)
  {
    art::flattenCollections(in, out);
    return true;
  }

  bool Mu2eProductMixer::mixStrawDigiADCWaveforms(std::vector<StrawDigiADCWaveformCollection const*> const& in,
                     StrawDigiADCWaveformCollection& out,
                     art::PtrRemapper const& remap)
  {
    art::flattenCollections(in, out);
    return true;
  }

  bool Mu2eProductMixer::mixStrawDigiMCs(std::vector<StrawDigiMCCollection const*> const& in,
                     StrawDigiMCCollection& out,
                     art::PtrRemapper const& remap)
  {
    art::flattenCollections(in, out);
    return true;
  }

  bool Mu2eProductMixer::mixEventWindowMarkers(std::vector<EventWindowMarker const*> const& in,
                     EventWindowMarker& out,
                     art::PtrRemapper const& remap){
    // assert that only one EventWindowMarker be present
    if (in.size() != 1){
      std::string msg = "attempting to mix more than 1 EventWindowMarker: ";
      msg += std::to_string(in.size()) + " present";
      throw cet::exception("BADINPUT") << msg << std::endl;
    }
    out = *in[0];
    return true;
  }

  //----------------------------------------------------------------
  bool Mu2eProductMixer::mixEventIDs(std::vector<art::EventIDSequence const*> const &in,
                                     art::EventIDSequence& out,
                                     art::PtrRemapper const&)
  {
    art::flattenCollections(in, out);
    return true;
  }

  //----------------------------------------------------------------
  bool Mu2eProductMixer::mixVolumeInfos(std::vector<PhysicalVolumeInfoMultiCollection const*> const &in,
                                        PhysicalVolumeInfoMultiCollection& out,
                                        art::PtrRemapper const&)
  {
    if(!in.empty()) {
      // We add incoming data to the smaller event-level structure that eliminates
      // some duplicates.  Then we transfer unuque event level data into the larger
      // subrun structure, again eliminating duplicates.

      const auto numStages = in[0]->size();
      std::vector<VolumeMap> eventInfos(numStages);

      for(const auto& mcoll: in) {

        if(mcoll->size() != numStages) {
          throw cet::exception("BADINPUT")<<"Mu2eProductMixer/evt: incompatible PhysicalVolumeInfoMultiCollection inputs. "
                                          <<"numStages="<<numStages<<" vs "<<mcoll->size()
                                          <<std::endl;
        }

        for(unsigned stage=0; stage<numStages; ++stage) {
          for(const auto& entry: (*mcoll)[stage]) {
            addInfo(&eventInfos[stage], entry);
          }
        }
      }

      out.clear();
      out.resize(numStages);
      for(unsigned stage=0; stage<numStages; ++stage) {
        out[stage].insert(eventInfos[stage].begin(), eventInfos[stage].end());
      }


      if(subrunVolumes_.empty()) {
        subrunVolumes_.resize(numStages);
      }
      else {
        if(subrunVolumes_.size() != numStages) {
          throw cet::exception("BADINPUT")<<"Mu2eProductMixer/sr: incompatible PhysicalVolumeInfoMultiCollection inputs. "
                                          <<"numStages="<<numStages<<" vs "<<subrunVolumes_.size()
                                          <<std::endl;
        }
      }

      for(unsigned stage=0; stage<numStages; ++stage) {
        for(const auto& entry: eventInfos[stage]) {
          addInfo(&subrunVolumes_[stage], entry);
        }
      }

    }

    const bool putVolsIntoEvent{evtVolInstanceName_};
    return putVolsIntoEvent;
  }

    //----------------------------------------------------------------
  bool Mu2eProductMixer::mixGenEventCount(std::vector<GenEventCount const*> const& in,
                                       GenEventCount& out,
                                       art::PtrRemapper const& remap)
  {
    if(in.size() > 1) {
        throw cet::exception("BADINPUT")<<"Mu2eProductMixer/subrun: can't mix GenEventCount" << std::endl;
    } else if(in.size() == 1) {
      generatedEvents_ = in[0]->count();
    }

    return false;
  }


  //----------------------------------------------------------------
  bool Mu2eProductMixer::mixCosmicLivetime(std::vector<CosmicLivetime const*> const& in,
                                                 CosmicLivetime& out,
                                                 art::PtrRemapper const& remap)
  {
    if(in.size() > 1) {
        throw cet::exception("BADINPUT")<<"Mu2eProductMixer/subrun: can't mix CosmicLiveTime" << std::endl;
    } else if(in.size() == 1) {
      area_ = in[0]->area();
      lowE_ = in[0]->lowE();
      highE_ = in[0]->highE();
      fluxConstant_ = in[0]->fluxConstant();
      totalPrimaries_ = in[0]->primaries();
      livetime_ = in[0]->liveTime();
    }

    return true;
  }

  //----------------------------------------------------------------
  void Mu2eProductMixer::addInfo(VolumeMap* map, const PhysicalVolumeInfoSingleStage::value_type& entry) {
    const auto it = map->find(entry.first);
    if(it != map->end()) {
      if(it->second != entry.second) {
        throw cet::exception("BADINPUT")<<"Mu2eProductMixer::addInfo(): inconsistent volume infos for index "
                                        <<it->first.asInt()<<": "
                                        <<"a = "<<it->second
                                        <<", b = "<<entry.second
                                        <<std::endl;
      }
    }
    else {
      map->insert(entry);
    }
  }

  //----------------------------------------------------------------

}
//================================================================
