// Andrei Gaponenko, 2018

#include "EventMixing/inc/Mu2eProductMixer.hh"

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
  Mu2eProductMixer::Mu2eProductMixer(const Config& conf, art::MixHelper& helper) {

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

    for(const auto& e: conf.caloShowerStepMixer().mixingMap()) {
      helper.declareMixOp
        (e.inTag, e.resolvedInstanceName(), &Mu2eProductMixer::mixCaloShowerSteps, *this);
    }

    for(const auto& e: conf.extMonSimHitMixer().mixingMap()) {
      helper.declareMixOp
        (e.inTag, e.resolvedInstanceName(), &Mu2eProductMixer::mixExtMonSimHits, *this);
    }

    for(const auto& e: conf.protonBunchIntensityMixer().mixingMap()) {
      helper.declareMixOp
        (e.inTag, e.resolvedInstanceName(), &Mu2eProductMixer::mixProtonBunchIntensity, *this);
    }

    for(const auto& e: conf.protonTimeMapMixer().mixingMap()) {
      helper.declareMixOp
        (e.inTag, e.resolvedInstanceName(), &Mu2eProductMixer::mixProtonTimeMap, *this);
    }

  }

  //----------------------------------------------------------------
  bool Mu2eProductMixer::mixGenParticles(std::vector<GenParticleCollection const*> const& in,
                                         GenParticleCollection& out,
                                         art::PtrRemapper const& remap)
  {
    art::flattenCollections(in, out, genOffsets_);
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
  bool Mu2eProductMixer::mixProtonBunchIntensity(std::vector<ProtonBunchIntensity const*> const& in,
                                                 ProtonBunchIntensity& out,
                                                 art::PtrRemapper const& remap)
  {
    for(const auto& x: in) {
      out.add(*x);
    }

    return true;
  }

  bool Mu2eProductMixer::mixProtonTimeMap(std::vector<SimParticleTimeMap const*> const& in,
                                                 SimParticleTimeMap& out,
                                                 art::PtrRemapper const& remap)
  {
   for(size_t incount = 0; incount < in.size(); ++incount) {
      auto const& timemap = *in[incount]; 
      //std::cout << "Mixing time map " << incount << " size " << timemap.size() << std::endl;
      for(auto & imap : timemap) {
	auto newptr = remap(imap.first, simOffsets_[incount]);
	out[newptr] = imap.second;
	// do I need to go down the chain?  I think not
      }
    }

    return true;
  }

}
//================================================================
