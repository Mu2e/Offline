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
#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "Mu2eUtilities/inc/SimParticleParentGetter.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoMultiCollection.hh"
#include "Mu2eUtilities/inc/PhysicalVolumeMultiHelper.hh"

namespace mu2e {

  //================================================================
  class CompressPhysicalVolumes : public art::EDProducer {
    art::InputTag volumesInput_;

    typedef std::vector<art::InputTag> InputTags;
    InputTags hitInputs_;
    InputTags particleInputs_;

    typedef PhysicalVolumeInfoSingleStage::key_type key_type;
    typedef std::vector<std::set<key_type> > UsedKeys;
    UsedKeys used_;

    const PhysicalVolumeInfoMultiCollection *incoll_;
    std::unique_ptr<PhysicalVolumeMultiHelper> helper_;

  public:
    explicit CompressPhysicalVolumes(const fhicl::ParameterSet& pset);
    virtual void produce(art::Event& event) override;
    virtual void beginSubRun(art::SubRun& sr) override;
    virtual void endSubRun(art::SubRun& sr) override;
  };

  //================================================================
  CompressPhysicalVolumes::CompressPhysicalVolumes(const fhicl::ParameterSet& pset)
    : volumesInput_(pset.get<std::string>("volumesInput"))
    , incoll_(nullptr)
    , helper_(nullptr)
  {
    produces<PhysicalVolumeInfoMultiCollection,art::InSubRun>();

    typedef std::vector<std::string> VS;
    const VS hi(pset.get<VS>("hitInputs"));
    for(const auto& i : hi) {
      hitInputs_.emplace_back(i);
    }

    const VS pi(pset.get<VS>("particleInputs"));
    for(const auto& i : pi) {
      particleInputs_.emplace_back(i);
    }
  }

  //================================================================
  void CompressPhysicalVolumes::beginSubRun(art::SubRun& sr) {
    art::Handle<PhysicalVolumeInfoMultiCollection> ih;
    sr.getByLabel(volumesInput_, ih);
    incoll_ = &*ih;
    helper_.reset(new PhysicalVolumeMultiHelper(*incoll_));

    used_.clear();
    used_.resize(incoll_->size());
  }

  //================================================================
  void CompressPhysicalVolumes::produce(art::Event& event) {

    // Go through all SimParticles in the specified products and
    // record their begin and end volume indexes.

    for(const auto& tag : hitInputs_) {
      auto ih = event.getValidHandle<StepPointMCCollection>(tag);
      for(const auto& hit : *ih) {
        const SimParticle& p = *hit.simParticle();
        const PhysicalVolumeInfoMultiCollection::size_type stage = helper_->iSimStage(p);
        used_[stage].insert(key_type(hit.simParticle()->startVolumeIndex()));
        used_[stage].insert(key_type(hit.simParticle()->endVolumeIndex()));
      }
    }

    for(const auto& tag : particleInputs_) {
      auto ih = event.getValidHandle<SimParticleCollection>(tag);
      for(const auto& spe : *ih) {
        const SimParticle& p = spe.second;
        const PhysicalVolumeInfoMultiCollection::size_type stage = helper_->iSimStage(p);
        used_[stage].insert(key_type(p.startVolumeIndex()));
        used_[stage].insert(key_type(p.endVolumeIndex()));
      }
    }
  }

  //================================================================
  void CompressPhysicalVolumes::endSubRun(art::SubRun& sr) {
    std::unique_ptr<PhysicalVolumeInfoMultiCollection> out(new PhysicalVolumeInfoMultiCollection);
    out->resize(incoll_->size());

    unsigned totalCount(0), passedCount(0);
    for(unsigned stage = 0; stage < incoll_->size(); ++stage) {
      (*out)[stage] = std::make_pair((*incoll_)[stage].first, PhysicalVolumeInfoSingleStage());
      const auto& ss = (*incoll_)[stage].second;
      totalCount += ss.size();
      for(const auto& in : ss) {
        if(used_[stage].find(in.first) != used_[stage].end()) {
          ++passedCount;
          (*out)[stage].second[in.first] = in.second;
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
