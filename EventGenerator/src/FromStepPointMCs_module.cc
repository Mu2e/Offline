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
#include "art/Persistency/Common/Assns.h"

#include "ConditionsService/inc/GlobalConstantsHandle.hh"
#include "ConditionsService/inc/ParticleDataTable.hh"
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
  }

  class FromStepPointMCs : public art::EDProducer {
  public:
    explicit FromStepPointMCs(fhicl::ParameterSet const& pset);
    void produce(art::Event& event);

  private:
    std::string inModuleLabel_;
    std::string  inInstanceName_;
    GlobalConstantsHandle<ParticleDataTable> pdt_;
    int logLevel_;
    bool allowDuplicates_;  // easily get a wrong answer by double counting particles
  };

  FromStepPointMCs::FromStepPointMCs(fhicl::ParameterSet const& pset)
    : inModuleLabel_(pset.get<std::string>("inputModuleLabel"))
    , inInstanceName_(pset.get<std::string>("inputInstanceName"))
    , logLevel_(pset.get<int>("logLevel", 0))
    , allowDuplicates_(pset.get<bool>("allowDuplicates", false))
  {
    produces<GenParticleCollection>();
    produces<GenParticleSPMHistory>();
  }

  void FromStepPointMCs::produce(art::Event& event) {

    std::auto_ptr<GenParticleCollection> output(new GenParticleCollection);
    std::auto_ptr<GenParticleSPMHistory> history(new GenParticleSPMHistory);

    art::ProductID gpc_pid = getProductID<GenParticleCollection>(event);

    // The input collection
    art::Handle<mu2e::StepPointMCCollection> ih;
    event.getByLabel(inModuleLabel_, inInstanceName_, ih);
    const StepPointMCCollection& inhits(*ih);

    // The purpose of this module is to continue simulation of an
    // event from saved hits (normally in Virtual Detectors).  If our
    // inputs contain multiple hits made by the same particle
    // then creation of new GenParticle from each hit would be wrong.
    // Make sure each SimParticle occurs only once.

    typedef std::set<const SimParticle*> SPSet;
    SPSet seenParticles;
    for(StepPointMCCollection::const_iterator i=inhits.begin(); i!=inhits.end(); ++i) {
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
                                      i->time()
                                      ));

        history->addSingle(art::Ptr<GenParticle>(gpc_pid, output->size()-1, event.productGetter(gpc_pid)),
                           art::Ptr<StepPointMC>(ih, std::distance(inhits.begin(), i))
                           );
      }
      else {
        mf::LogWarning warn("FromStepPointMCs");
        warn<<"WARNING: No particle data for pdgId = "<<particle->pdgId()<<"\n";
      }

    }

    event.put(std::move(output));
    event.put(std::move(history));
  }
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::FromStepPointMCs);
