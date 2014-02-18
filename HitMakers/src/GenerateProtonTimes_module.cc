// Draw time values from the proton time profile distribution and
// associate SimParticles with time offsets of their primary protons.
//
// Andrei Gaponenko, 2013

#include <string>
#include <memory>

#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "CLHEP/Random/RandGaussQ.h"

#include "art/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Utilities/InputTag.h"

#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/SimParticleTimeMap.hh"
#include "SeedService/inc/SeedService.hh"
#include "Mu2eUtilities/inc/ProtonPulseRandPDF.hh"
#include "Mu2eUtilities/inc/SimParticleCollectionPrinter.hh"

namespace mu2e {

  class GenerateProtonTimes : public art::EDProducer {

  public:

    explicit GenerateProtonTimes(const fhicl::ParameterSet& pset);

    virtual void produce(art::Event& e) override;

  private:
    ProtonPulseRandPDF  protonPulse_;
    int  verbosityLevel_;
  };

  //================================================================
  GenerateProtonTimes::GenerateProtonTimes(fhicl::ParameterSet const& pset)
    : protonPulse_(createEngine(art::ServiceHandle<SeedService>()->getSeed()), 
                   pset.get<std::string>("pulseType","default") )
    , verbosityLevel_(pset.get<int>("verbosityLevel", 0))
  {
    produces<SimParticleTimeMap>();
  }

  //================================================================
  void GenerateProtonTimes::produce(art::Event& event) {
    std::unique_ptr<SimParticleTimeMap> res(new SimParticleTimeMap);

    std::vector<art::Handle<SimParticleCollection> > colls;
    event.getManyByType(colls);

    // Generate and record offsets for all primaries
    for(const auto& ih : colls) {
      for(const auto& iter : *ih) {
        if(iter.second.isPrimary()) {
          art::Ptr<SimParticle> part(ih, iter.first.asUint());
          if(!part->genParticle()->generatorId().isCosmic()) {
            (*res)[part] = protonPulse_.fire();
          }
          else {
            (*res)[part] = 0;
          }
        }
      }
    }

    if(verbosityLevel_ > 10) {
      std::cout<<"GenerateProtonTimes dump begin"<<std::endl;
      for(const auto& i : *res) {
        SimParticleCollectionPrinter::print(std::cout, *i.first);
        std::cout<<" => "<<i.second<<std::endl;
      }
      std::cout<<"GenerateProtonTimes dump end"<<std::endl;
    }

    event.put(std::move(res));
  } // end GenerateProtonTimes::produce.

  //================================================================
} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::GenerateProtonTimes)
