// Compresses  PhysicalVolumeInfoMultiCollection in SubRun
// by eliminating entries not referenced by any of SimParticles
// in a specified set of products.
//
// Andrei Gaponenko, 2013

#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <memory>

// art includes.
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "Mu2eUtilities/inc/SimParticleParentGetter.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoMultiCollection.hh"

namespace mu2e {

  //================================================================
  class CompressPhysicalVolumes : public art::EDProducer {
    art::ProductToken<PhysicalVolumeInfoMultiCollection> const volumesToken_;
    std::vector<art::ProductToken<StepPointMCCollection>> hitTokens_;
    std::vector<art::ProductToken<SimParticleCollection>> particleTokens_;

    typedef PhysicalVolumeInfoSingleStage::key_type key_type;
    typedef std::vector<std::set<key_type> > UsedKeys;
    UsedKeys used_;

    PhysicalVolumeInfoMultiCollection const* incoll_{nullptr};

    void produce(art::Event& event) override;
    void beginSubRun(art::SubRun& sr) override;
    void endSubRun(art::SubRun& sr) override;

  public:

    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::Atom<art::InputTag> volumesInput {
        Name("volumesInput"),
          Comment("The collection to compress")
          };

      fhicl::Sequence<art::InputTag> hitInputs {
        Name("hitInputs"),
          Comment("A list of StepPointMCCollections; physical volumes referred to by corresponding SimParticles will be kept.")
          };

      fhicl::Sequence<art::InputTag> particleInputs {
        Name("particleInputs"),
          Comment("A list of SimParticleCollections; physical volumes referred to by the particles will be kept.")
          };
    };

    using Parameters = art::EDProducer::Table<Config>;
    explicit CompressPhysicalVolumes(const Parameters& conf);
  };

  //================================================================
  CompressPhysicalVolumes::CompressPhysicalVolumes(const Parameters& conf)
    : art::EDProducer{conf},
    volumesToken_{consumes<PhysicalVolumeInfoMultiCollection, art::InSubRun>(conf().volumesInput())}
  {
    produces<PhysicalVolumeInfoMultiCollection,art::InSubRun>();

    for(const auto& i : conf().hitInputs()) {
      hitTokens_.emplace_back(consumes<StepPointMCCollection>(i));
    }

    for(const auto& i : conf().particleInputs()) {
      particleTokens_.emplace_back(consumes<SimParticleCollection>(i));
    }
  }

  //================================================================
  void CompressPhysicalVolumes::beginSubRun(art::SubRun& sr) {
    auto const& ih = sr.getValidHandle(volumesToken_);
    incoll_ = &*ih;
    used_.clear();
    used_.resize(incoll_->size());
  }

  //================================================================
  void CompressPhysicalVolumes::produce(art::Event& event)
  {
    // Go through all SimParticles in the specified products and
    // record their begin and end volume indexes.

    for(const auto& token : hitTokens_) {
      auto ih = event.getValidHandle(token);
      for(const auto& hit : *ih) {
        const SimParticle& p = *hit.simParticle();
        const PhysicalVolumeInfoMultiCollection::size_type stage = p.simStage();
        used_[stage].insert(key_type(hit.simParticle()->startVolumeIndex()));
        used_[stage].insert(key_type(hit.simParticle()->endVolumeIndex()));
      }
    }

    for(const auto& token : particleTokens_) {
      auto ih = event.getValidHandle(token);
      for(const auto& spe : *ih) {
        const SimParticle& p = spe.second;
        const PhysicalVolumeInfoMultiCollection::size_type stage = p.simStage();
        used_[stage].insert(key_type(p.startVolumeIndex()));
        used_[stage].insert(key_type(p.endVolumeIndex()));
      }
    }
  }

  //================================================================
  void CompressPhysicalVolumes::endSubRun(art::SubRun& sr)
  {
    auto out = std::make_unique<PhysicalVolumeInfoMultiCollection>();
    out->resize(incoll_->size());

    unsigned totalCount(0), passedCount(0);
    for(unsigned stage = 0; stage < incoll_->size(); ++stage) {
      (*out)[stage] = PhysicalVolumeInfoSingleStage();
      const auto& ss = (*incoll_)[stage];
      totalCount += ss.size();
      for(const auto& in : ss) {
        if(used_[stage].find(in.first) != used_[stage].end()) {
          ++passedCount;
          (*out)[stage][in.first] = in.second;
        }
      }
    }

    sr.put(std::move(out));

    mf::LogInfo("Summary")
      << "CompressPhysicalVolumes stats: passed "
      << passedCount <<" / "<<totalCount<<" volume entries\n";
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::CompressPhysicalVolumes);
