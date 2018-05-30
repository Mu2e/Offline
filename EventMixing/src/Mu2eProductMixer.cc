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
    std::cout<<"Mu2eProductMixer ctr begin"<<std::endl;

    for(const auto& e: conf.genParticleMixer().mixingMap()) {
      std::cout<<"will mix GenParticles "<<e.inTag<<"  =====>  "<<e.outInstance<<std::endl;
      helper.declareMixOp
        (e.inTag, e.outInstance, &Mu2eProductMixer::mixGenParticles, *this);
    }

    for(const auto& e: conf.simParticleMixer().mixingMap()) {
      std::cout<<"will mix SimParticles "<<e.inTag<<"  =====>  "<<e.outInstance<<std::endl;
      helper.declareMixOp
        (e.inTag, e.outInstance, &Mu2eProductMixer::mixSimParticles, *this);
    }

    for(const auto& e: conf.stepPointMCMixer().mixingMap()) {
      std::cout<<"will mix StepPointMCs "<<e.inTag<<"  =====>  "<<e.outInstance<<std::endl;
      helper.declareMixOp
        (e.inTag, e.outInstance, &Mu2eProductMixer::mixStepPointMCs, *this);
    }

    for(const auto& e: conf.caloShowerStepMixer().mixingMap()) {
      std::cout<<"will mix CaloShowerSteps "<<e.inTag<<"  =====>  "<<e.outInstance<<std::endl;
      helper.declareMixOp
        (e.inTag, e.outInstance, &Mu2eProductMixer::mixCaloShowerSteps, *this);
    }

    for(const auto& e: conf.extMonSimHitMixer().mixingMap()) {
      std::cout<<"will mix ExtMonSimHits "<<e.inTag<<"  =====>  "<<e.outInstance<<std::endl;
      helper.declareMixOp
        (e.inTag, e.outInstance, &Mu2eProductMixer::mixExtMonSimHits, *this);
    }

    std::cout<<"Mu2eProductMixer ctr end"<<std::endl;
  }

  //----------------------------------------------------------------
  bool Mu2eProductMixer::mixGenParticles(std::vector<GenParticleCollection const*> const& in,
                                         GenParticleCollection& out,
                                         art::PtrRemapper const& remap)
  {
    std::cout<<"Mu2eProductMixer::mixGenParticles() BEGIN\n";

    art::flattenCollections(in, out, genOffsets_);

    std::cout<<"genOffsets = ";
    std::copy(genOffsets_.begin(), genOffsets_.end(), std::ostream_iterator<unsigned>(std::cout, " "));
    std::cout<<std::endl;

    std::cout<<"Mu2eProductMixer::mixGenParticles() END\n";
    return true;
  }

  //----------------------------------------------------------------
  bool Mu2eProductMixer::mixSimParticles(std::vector<SimParticleCollection const*> const& in,
                                         SimParticleCollection& out,
                                         art::PtrRemapper const& remap)
  {
    // Flatten the input collections; does not update Ptrs.
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

}
//================================================================
