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
#include "MCDataProducts/inc/SimParticleProtonPulseTimeMap.hh"
#include "SeedService/inc/SeedService.hh"
#include "Mu2eUtilities/inc/ProtonPulseRandPDF.hh"
#include "Mu2eUtilities/inc/SimParticleCollectionPrinter.hh"

namespace mu2e {

  namespace {
    // In art v1_00_06 the ValidHandle::id() method needed by an
    // art::Ptr constructor is missing.  A workaround:
    template<class PROD> struct MyHandle : public art::ValidHandle<PROD> {
      MyHandle(const art::ValidHandle<PROD>& h) : art::ValidHandle<PROD>(h) {}
      art::ProductID id( ) const { return this->provenance()->productID(); }
    };
    template<class PROD> MyHandle<PROD> makeMyHandle(const art::ValidHandle<PROD>& h) {
      return MyHandle<PROD>(h);
    }

  } // end of anonymous namespace

  class GenerateProtonTimes : public art::EDProducer {

  public:

    explicit GenerateProtonTimes(const fhicl::ParameterSet& pset);

    virtual void produce(art::Event& e) override;

  private:
    ProtonPulseRandPDF  protonPulse_;
    art::InputTag simParticlesTag_;
    bool mapDaughters_;
    int  verbosityLevel_;
  };

  //================================================================
  GenerateProtonTimes::GenerateProtonTimes(fhicl::ParameterSet const& pset)
    : protonPulse_(createEngine(art::ServiceHandle<SeedService>()->getSeed()))
    , simParticlesTag_(pset.get<std::string>("simParticlesTag"))
    , mapDaughters_(pset.get<bool>("mapDaughters", false))
    , verbosityLevel_(pset.get<int>("verbosityLevel", 0))
  {
    produces<SimParticleProtonPulseTimeMap>();
  }

  //================================================================
  void GenerateProtonTimes::produce(art::Event& event) {
    std::unique_ptr<SimParticleProtonPulseTimeMap> res(new SimParticleProtonPulseTimeMap);

    auto ih = event.getValidHandle<SimParticleCollection>(simParticlesTag_);

    // Generate and record offsets for all primaries
    for(const auto& iter : *ih) {
      if(iter.second.isPrimary()) {
        art::Ptr<SimParticle> part(makeMyHandle(ih), iter.first.asUint());
        (*res)[part] = protonPulse_.fire();
      }
    }

    if(mapDaughters_) {
      // Associate the offsets to the daughters as well
      for(const auto& iter : *ih) {
        if(!iter.second.isPrimary()) {
          art::Ptr<SimParticle> daughter(makeMyHandle(ih), iter.first.asUint());

          art::Ptr<SimParticle> mother = daughter->parent();
          while(mother && !mother->isPrimary()) {
            mother = mother->parent();
          }

          const auto pm = res->find(mother);
          if(pm == res->end()) {
            throw cet::exception("BADINPUTS")<<"ERROR: GenerateProtonTimes: a particle primary is not available\n";
          }

          (*res)[daughter] = pm->second;
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
