// Identify stopped particles in a given input collection
// and write them to a new SimParticlePtrCollection.
//
// Andrei Gaponenko, 2013

#include <string>
#include <vector>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <limits>

#include "cetlib_except/exception.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"

// Mu2e includes.
#include "DataProducts/inc/PDGCode.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/SimParticlePtrCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoMultiCollection.hh"
#include "Mu2eUtilities/inc/PhysicalVolumeMultiHelper.hh"

#include "TH1D.h"

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
    std::vector<std::string> vetoedMaterials_;

    int verbosityLevel_;

    unsigned offsetThreshold_; // to select particles from the current simulation stage

    typedef std::set<PDGCode::type> PDGCodeSet;
    PDGCodeSet particleTypes_;

    TH1* hStopMaterials_;

    const PhysicalVolumeInfoMultiCollection *vols_;

    bool isStopped(const SimParticle& particle) const;
    bool materialAccepted(const std::string& material) const;

    unsigned numTotalParticles_;
    unsigned numRequestedTypeStops_;
    unsigned numRequestedMateralStops_;
  };

  //================================================================
  StoppedParticlesFinder::StoppedParticlesFinder(const fhicl::ParameterSet& pset)
    : particleInput_(pset.get<std::string>("particleInput"))
    , physVolInfoInput_(pset.get<std::string>("physVolInfoInput"))
    , stoppingMaterial_(pset.get<std::string>("stoppingMaterial", ""))
    , vetoedMaterials_(pset.get<std::vector<std::string> >("vetoedMaterials", std::vector<std::string>()))
    , verbosityLevel_(pset.get<int>("verbosityLevel", 0))
    , offsetThreshold_(pset.get<unsigned>("particleOffsetThreshold", 0))
    , hStopMaterials_(art::ServiceHandle<art::TFileService>()->make<TH1D>("stopmat", "Stopping materials", 1, 0., 1.))
    , vols_()
    , numTotalParticles_()
    , numRequestedTypeStops_()
    , numRequestedMateralStops_()
  {
    produces<SimParticlePtrCollection>();

    auto pt(pset.get<std::vector<int> >("particleTypes"));
    for(const auto& pid : pt) {
      particleTypes_.insert(PDGCode::type(pid));
    }

    //----------------
    if(verbosityLevel_ > 0) {
      std::ostringstream os;
      os<<"Particle types: [ ";
      std::copy(particleTypes_.begin(), particleTypes_.end(), std::ostream_iterator<int>(os, ", "));
      os<<" ]"<<std::endl;

      os<<"offsetThreshold  = "<<offsetThreshold_<<std::endl;

      os<<"stoppingMaterial = "<<stoppingMaterial_<<std::endl;

      os<<"vetoedMaterials = [ ";
      std::copy(vetoedMaterials_.begin(), vetoedMaterials_.end(), std::ostream_iterator<std::string>(os, ", "));
      os<<" ]"<<std::endl;

      mf::LogInfo("Info")<<os.str();
    }

    //----------------
    if( !stoppingMaterial_.empty() &&  !vetoedMaterials_.empty()) {
      throw cet::exception("BADCONFIG")<<"Only one stoppingMaterial or vetoedMaterials may be requested.";
    }

  }

  //================================================================
  void StoppedParticlesFinder::beginSubRun(art::SubRun& sr) {
    art::Handle<PhysicalVolumeInfoMultiCollection> volh;
    sr.getByLabel(physVolInfoInput_, volh);
    vols_ = &*volh;

    if(verbosityLevel_ > 1) {
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
      if((particle.id().asUint() >= offsetThreshold_)
         &&(particleTypes_.find(particle.pdgId()) != particleTypes_.end())
         && isStopped(particle))
        {
          ++numRequestedTypeStops_;

          if(verbosityLevel_ > 2) {
            std::cout<<"stopped particle "<<particle.pdgId()
                     <<" at pos="<<particle.endPosition()
                     <<" time="<<particle.endGlobalTime()
                     <<" reason="<<particle.stoppingCode()
                     <<" in volume "<<vi.endVolume(particle)
                     <<std::endl;
          }

          const std::string material = vi.endVolume(particle).materialName();
          hStopMaterials_->Fill(material.c_str(), 1.);

          if(materialAccepted(material)) {
            ++numRequestedMateralStops_;
            output->emplace_back(ih, particle.id().asUint());
          }
        }
    }

    event.put(std::move(output));
  }

  //================================================================
  bool StoppedParticlesFinder::isStopped(const SimParticle& particle) const {
    return particle.endMomentum().v().mag2() <= std::numeric_limits<double>::epsilon();
  }

  //================================================================
  bool StoppedParticlesFinder::materialAccepted(const std::string& material) const {
    // Check if the stop is in a material of interest.  The string
    // comparison looks inefficient, however this is called only for
    // stopped particles that are a small fraction of all particles.
    // One could pre-compute a lookup table for stopping volume index,
    // (rather than the material name string) but that would require
    // going through all of the O(10^4) volumes at the initialization
    // stage.

    bool ret = false;
    if(stoppingMaterial_.empty()) {
      ret = true;
      for(const auto& i : vetoedMaterials_) {
        if(i == material) {
          ret = false;
        }
      }
    }
    else {
      ret = (material == stoppingMaterial_);
    }

    return ret;
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
