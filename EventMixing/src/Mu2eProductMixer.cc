// Andrei Gaponenko, 2018

#include "Offline/EventMixing/inc/Mu2eProductMixer.hh"

#include <utility>
#include <algorithm>
#include <iterator>
#include <iostream>

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
      , applyTimeOffset_{! conf.simTimeOffset().empty() }
      , timeOffsetTag_{ conf.simTimeOffset() }
      , stoff_(0.0)
  {
    if(applyTimeOffset_){
      std::cout << "Mu2eProductMixer: Applying time offsets from " << timeOffsetTag_ << std::endl;
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

    for(const auto& e: conf.cosmicLivetimeMixer().mixingMap()) {
      helper.declareMixOp
        (e.inTag, e.resolvedInstanceName(), &Mu2eProductMixer::mixCosmicLivetime, *this);
    }

    for(const auto& e: conf.eventIDMixer().mixingMap()) {
      helper.declareMixOp
        (e.inTag, e.resolvedInstanceName(), &Mu2eProductMixer::mixEventIDs, *this);
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
  }

  void Mu2eProductMixer::startEvent(art::Event const& e) {
    if(applyTimeOffset_){
    // find the time offset in the event, and copy it locally
      const auto& stoH = e.getValidHandle<SimTimeOffset>(timeOffsetTag_);
      stoff_ = *stoH;
    }
  }


  //================================================================
  void Mu2eProductMixer::beginSubRun(const art::SubRun&) {
    subrunVolumes_.clear();
  }

  //----------------------------------------------------------------
  void Mu2eProductMixer::endSubRun(art::SubRun& sr) {
    if(mixVolumes_) {
      auto col =  std::make_unique<PhysicalVolumeInfoMultiCollection>();

      col->resize(subrunVolumes_.size());
      for(unsigned stage=0; stage<subrunVolumes_.size(); ++stage) {
        (*col)[stage].insert(subrunVolumes_[stage].begin(), subrunVolumes_[stage].end());
      }

      sr.put(std::move(col), subrunVolInstanceName_);
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

  //----------------------------------------------------------------
  bool Mu2eProductMixer::mixCosmicLivetime(std::vector<CosmicLivetime const*> const& in,
                                                 CosmicLivetime& out,
                                                 art::PtrRemapper const& remap)
  {
    if(in.size() > 1)
      throw cet::exception("BADINPUT")<<"Mu2eProductMixer/evt: can't mix CosmicLiveTime" << std::endl; 
    for(const auto& x: in) {
      out = *x;
    }

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
