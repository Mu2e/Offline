// Identify stopped particles in a given input collection
// and write them to a new SimParticlePtrCollection.
//
// Andrei Gaponenko, 2013

#include <string>
#include <vector>
#include <limits>

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"

// Mu2e includes.
#include "MCDataProducts/inc/PDGCode.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/SimParticlePtrCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoMultiCollection.hh"
#include "Mu2eUtilities/inc/PhysicalVolumeMultiHelper.hh"

namespace mu2e {

  class StoppedParticlesFinder : public art::EDProducer {
  public:
    explicit StoppedParticlesFinder(fhicl::ParameterSet const& pset);
    void beginSubRun(art::SubRun& sr) override;
    void produce(art::Event& evt) override;
    void endJob() override;
  private:
    art::InputTag particleInput_;
    art::InputTag physVolInfoInput_;
    std::string stoppingMaterial_;
    int verbosityLevel_;

    typedef std::set<PDGCode::type> PDGCodeSet;
    PDGCodeSet particleTypes_;

    const PhysicalVolumeInfoMultiCollection *vols_;

    bool isStopped(const SimParticle& particle);

    unsigned numTotalParticles_;
    unsigned numRequestedTypeStops_;
    unsigned numRequestedMateralStops_;
  };

  //================================================================
  StoppedParticlesFinder::StoppedParticlesFinder(const fhicl::ParameterSet& pset)
    : particleInput_(pset.get<std::string>("particleInput"))
    , physVolInfoInput_(pset.get<std::string>("physVolInfoInput"))
    , stoppingMaterial_(pset.get<std::string>("stoppingMaterial"))
    , verbosityLevel_(pset.get<int>("verbosityLevel", 0))
    , vols_()
    , numTotalParticles_()
    , numRequestedTypeStops_()
    , numRequestedMateralStops_()
  {
    auto pt(pset.get<std::vector<int> >("particleTypes"));
    for(const auto& pid : pt) {
      particleTypes_.insert(PDGCode::type(pid));
    }

    produces<SimParticlePtrCollection>();
  }

  //================================================================
  void StoppedParticlesFinder::beginSubRun(art::SubRun& sr) {
    art::Handle<PhysicalVolumeInfoMultiCollection> volh;
    sr.getByLabel(physVolInfoInput_, volh);
    vols_ = &*volh;

    if(verbosityLevel_ > 0) {
      std::cout<<"PhysicalVolumeInfoMultiCollection dump begin"<<std::endl;
      for(const auto& i : *vols_) {
        std::cout<<"*********************************************************"<<std::endl;
        std::cout<<"SimParticleNumberOffset = "<<i.first<<", collection size = "<<i.second.size()<<std::endl;
        for(const auto& entry : i.second) {
          std::cout<<entry.second<<std::endl;
        }
      }
      std::cout<<"PhysicalVolumeInfoMultiCollection dump end"<<std::endl;
    }
  }

  //================================================================
  void StoppedParticlesFinder::produce(art::Event& event) {

    std::unique_ptr<SimParticlePtrCollection> output(new SimParticlePtrCollection());

    PhysicalVolumeMultiHelper vi(*vols_);
    auto ih = event.getValidHandle<SimParticleCollection>(particleInput_);
    numTotalParticles_ += ih->size();
    for(const auto& i : *ih) {
      const SimParticle& particle = i.second;
      if((particleTypes_.find(particle.pdgId()) != particleTypes_.end())
         && isStopped(particle))
        {
          ++numRequestedTypeStops_;

          if(verbosityLevel_ > 1) {
            std::cout<<"stopped particle "<<particle.pdgId()
                     <<" at pos="<<particle.endPosition()
                     <<" time="<<particle.endGlobalTime()
                     <<" reason="<<particle.stoppingCode()
                     <<" in volume "<<vi.endVolume(particle)
                     <<std::endl;
          }

          // Check if the stop is in a material of interest.  The
          // string comparison in this loop looks inefficient, however
          // only a small fraction of particles are stopped.  One
          // could pre-compute a lookup table for stopping volume
          // index, but that would require going through all of the
          // O(10^4) volumes at the initialization stage.
          if(stoppingMaterial_.empty() ||
             (vi.endVolume(particle).materialName() == stoppingMaterial_))
            {
              ++numRequestedMateralStops_;
              output->emplace_back(ih, particle.id().asUint());
            }
        }
    }

    event.put(std::move(output));
  }

  //================================================================
  bool StoppedParticlesFinder::isStopped(const SimParticle& particle) {
    return particle.endMomentum().v().mag2() <= std::numeric_limits<double>::epsilon();
  }

  //================================================================
  void StoppedParticlesFinder::endJob() {
    mf::LogInfo("Summary")
      <<"StoppedParticlesFinder stats:"
      <<" accepted = "<<numRequestedMateralStops_
      <<", requested type stops = "<<numRequestedTypeStops_
      <<", total input particles = "<<numTotalParticles_
      << "\n";
  }

  //================================================================

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::StoppedParticlesFinder);
