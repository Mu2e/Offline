// Read a StepPointMC collection and create a GenParticleCollection
// from the recorded hits.
//
// Andrei Gaponenko, 2012

#include <iostream>
#include <string>
#include <memory>

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/Assns.h"

#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/GenParticleSPMHistory.hh"

namespace mu2e {

  namespace {
    std::string toStr(const SimParticle& particle) {
      std::ostringstream os;
      os<<"SimParticle("
        <<"id = "<<particle.id()
        <<", pdgId="<<particle.pdgId()
        <<", hasParent="<<particle.hasParent()
        <<", parentId = "<<particle.parentId()
        <<", startPosition="<<particle.startPosition()
        <<", endPosition="<<particle.endPosition()
        <<")";
      return os.str();
    }


    // In art v1_00_06 the ValidHandle::id() method needed by an
    // art::Ptr constructor was missing.  A workaround:
    template<class PROD> struct MyHandle : public art::ValidHandle<PROD> {
      MyHandle(const art::ValidHandle<PROD>& h) : art::ValidHandle<PROD>(h) {}
      art::ProductID id( ) const { return this->provenance()->productID(); }
    };
    template<class PROD> MyHandle<PROD> makeMyHandle(const art::ValidHandle<PROD>& h) {
      return MyHandle<PROD>(h);
    }
  }

  class FromStepPointMCs : public art::EDProducer {
  public:
    explicit FromStepPointMCs(fhicl::ParameterSet const& pset);
    void produce(art::Event& event) override;

  private:
    typedef std::vector<art::InputTag> InputTags;
    InputTags inputs_;
    GlobalConstantsHandle<ParticleDataTable> pdt_;
    int logLevel_;
    bool allowDuplicates_;  // easily get a wrong answer by double counting particles
  };

  FromStepPointMCs::FromStepPointMCs(fhicl::ParameterSet const& pset)
    : art::EDProducer{pset}
    , logLevel_(pset.get<int>("logLevel", 0))
    , allowDuplicates_(pset.get<bool>("allowDuplicates", false))
  {
    produces<GenParticleCollection>();
    produces<GenParticleSPMHistory>();

    typedef std::vector<std::string> Strings;
    Strings tagstr(pset.get<Strings>("inputTags", Strings()));

    if(tagstr.empty()) { // legacy support
      mf::LogWarning warn("FromStepPointMCs");
      warn<<"WARNING: FromStepPointMCs: the inputTags list is not set or empty. "
          <<"Will try to use inputModuleLabel and inputInstanceName instead.\n";

      std::string label(pset.get<std::string>("inputModuleLabel"));
      std::string instance(pset.get<std::string>("inputInstanceName"));
      inputs_.emplace_back(label, instance);
    }
    else {
      for(const auto& s : tagstr) {
        inputs_.emplace_back(s);
      }
    }
  }

  void FromStepPointMCs::produce(art::Event& event) {

    std::unique_ptr<GenParticleCollection> output(new GenParticleCollection);
    std::unique_ptr<GenParticleSPMHistory> history(new GenParticleSPMHistory);

    art::ProductID gpc_pid = event.getProductID<GenParticleCollection>();

    for(const auto& tag : inputs_) {
      auto ih = event.getValidHandle<mu2e::StepPointMCCollection>(tag);

      // The purpose of this module is to continue simulation of an
      // event from saved hits (normally in Virtual Detectors).  If our
      // inputs contain multiple hits made by the same particle
      // then creation of new GenParticle from each hit would be wrong.
      // Make sure each SimParticle occurs only once.
      typedef std::set<const SimParticle*> SPSet;
      SPSet seenParticles;
      for(StepPointMCCollection::const_iterator i=ih->begin(); i!=ih->end(); ++i) {
        const art::Ptr<SimParticle>& particle = i->simParticle();
        if(!seenParticles.insert(&*particle).second) {
          std::ostringstream os;
          os <<"FromStepPointMCs: duplicate SimParticle "<<toStr(*particle)<<" in input hits collection. Hit: "<<*i<<" in event "<<event.id();
          if(allowDuplicates_) {
            mf::LogWarning warn("FromStepPointMCs");
            warn<<"WARNING: "<<os.str()<<"\n";
          }
          else {
            throw cet::exception("BADINPUTS")<<os.str();
          }
        }

        if(logLevel_ > 1) {
          std::cout<<"AG: creating new GenParticle from hit "<<*i<<std::endl;
        }

        if(pdt_->particle(particle->pdgId())) {
          const double mass = pdt_->particle(particle->pdgId()).ref().mass().value();
          const CLHEP::HepLorentzVector fourMom(i->momentum(), sqrt(mass*mass + i->momentum().mag2()));

          output->push_back(GenParticle(particle->pdgId(),
                                        GenId::fromStepPointMCs,
                                        i->position(),
                                        fourMom,
                                        i->time(),
                                        i->properTime()
                                        ));

          history->addSingle(art::Ptr<GenParticle>(gpc_pid, output->size()-1, event.productGetter(gpc_pid)),
                             art::Ptr<StepPointMC>(makeMyHandle(ih), std::distance(ih->begin(), i))
                             );
        }
        else {
          mf::LogWarning warn("FromStepPointMCs");
          warn<<"WARNING: No particle data for pdgId = "<<particle->pdgId()<<"\n";
        }

      } // for (StepPointMCCollection entries)
    } // for(inputs_)

    event.put(std::move(output));
    event.put(std::move(history));
  }
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::FromStepPointMCs);
