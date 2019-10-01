// Read a StepPointMC collection and create a GenParticleCollection
// from the recorded hits, rotate around target axis with a random
// phi angle
//
// Zhengyun You, 2013

#include <iostream>
#include <string>
#include <memory>
#include <vector>

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "canvas/Persistency/Common/Assns.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"

#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/GenParticleSPMHistory.hh"
#include "SeedService/inc/SeedService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "ProductionTargetGeom/inc/ProductionTarget.hh"

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

    typedef std::vector<int> Vpdg;
    typedef std::vector<std::string> Vstring;

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

  class FromStepPointMCsRotateTarget : public art::EDProducer {
  public:
    explicit FromStepPointMCsRotateTarget(fhicl::ParameterSet const& pset);
    void produce(art::Event& event) override;

  private:
    typedef std::vector<art::InputTag> InputTags;
    InputTags inputs_;
    Vpdg inPdgId_;
    GlobalConstantsHandle<ParticleDataTable> pdt_;
    int logLevel_;
    bool allowDuplicates_;  // easily get a wrong answer by double counting particles

    double phiMin_;   // rotate an angle between phiMin_ and phiMax_
    double phiMax_;

    art::RandomNumberGenerator::base_engine_t& eng_;
    CLHEP::RandFlat randFlat_;

    CLHEP::HepRotation targetRotation_; // rotates target frame to Mu2e frame
    CLHEP::Hep3Vector  targetCenter_;

    bool firstEvent_;
  };

  FromStepPointMCsRotateTarget::FromStepPointMCsRotateTarget(fhicl::ParameterSet const& pset)
    : EDProducer{pset}
    , inPdgId_(pset.get<Vpdg>("inputPdgIds", Vpdg()))
    , logLevel_(pset.get<int>("logLevel", 0))
    , allowDuplicates_(pset.get<bool>("allowDuplicates", false))
    , phiMin_(pset.get<double>("phiMin", 0.))
    , phiMax_(pset.get<double>("phiMax", 360.))
    , eng_(createEngine(art::ServiceHandle<SeedService>()->getSeed()))
    , randFlat_(eng_)
    , firstEvent_(true)
  {
    produces<GenParticleCollection>();
    produces<GenParticleSPMHistory>();

    typedef std::vector<std::string> Strings;
    Strings tagstr(pset.get<Strings>("inputTags", Strings()));

    if(tagstr.empty()) { // legacy support
      mf::LogWarning warn("FromStepPointMCsRotateTarget");
      warn<<"WARNING: FromStepPointMCsRotateTarget: the inputTags list is not set or empty. "
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

    phiMin_   *= CLHEP::deg;
    phiMax_   *= CLHEP::deg;

  }

  void FromStepPointMCsRotateTarget::produce(art::Event& event) {

    if (firstEvent_) {
      GeomHandle<ProductionTarget> target;
      targetRotation_ = target->protonBeamRotation();
      targetCenter_ = target->position();
      if(logLevel_ > 1) {
        std::cout << "targetCenter   " << targetCenter_ << std::endl;
        std::cout << "targetRotation " << targetRotation_ << std::endl;
      }
      firstEvent_ = false;
    }

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
        if(logLevel_ > 1) {
          std::cout << "particle id " << particle->pdgId() << " pos " << i->position() << std::endl;
        }
        Vpdg::iterator pdgFinder = find(inPdgId_.begin(),inPdgId_.end(), particle->pdgId());
        if (pdgFinder == inPdgId_.end()) {
          continue;
        }

        if(!seenParticles.insert(&*particle).second) {
          std::ostringstream os;
          os <<"FromStepPointMCsRotateTarget: duplicate SimParticle "<<toStr(*particle)<<" in input hits collection. Hit: "<<*i<<" in event "<<event.id();
          if(allowDuplicates_) {
            mf::LogWarning warn("FromStepPointMCsRotateTarget");
            warn<<"WARNING: "<<os.str()<<"\n";
          }
          else {
            continue;
            //throw cet::exception("BADINPUTS")<<os.str();
          }
        }

        if(logLevel_ > 1) {
          std::cout<<"creating new GenParticle from hit "<<*i<<std::endl;
        }

        if(pdt_->particle(particle->pdgId())) {
          const CLHEP::Hep3Vector pos = i->position();
          const CLHEP::Hep3Vector mom = i->momentum();
          if(logLevel_ > 1) std::cout << "pos " << pos << " mom " << mom << std::endl;

          CLHEP::Hep3Vector mom_tgt = targetRotation_.inverse()*mom;
          CLHEP::Hep3Vector pos_tgt = targetRotation_.inverse()*(pos-targetCenter_);
          if(logLevel_ > 1) std::cout << "pos_tgt " << pos_tgt << " mom_tgt " << mom_tgt << std::endl;

          // rotate by an random angle around phi and get new p direction in target frame
          double rotatePhi = randFlat_.fire()*(phiMax_-phiMin_);
          if(logLevel_ > 1) std::cout << "rotatePhi " << rotatePhi << std::endl;

          mom_tgt.rotateZ(rotatePhi);
          pos_tgt.rotateZ(rotatePhi);
          if(logLevel_ > 1) std::cout << "pos_tgt_rot " << pos_tgt << " mom_tgt_rot " << mom_tgt << std::endl;

          // rotate around y back to mu2e frame
          CLHEP::Hep3Vector mom_new = targetRotation_*mom_tgt;
          CLHEP::Hep3Vector pos_new = targetRotation_*pos_tgt+targetCenter_;
          if(logLevel_ > 1) std::cout << "pos_rot " << pos_new << " mom_rot " << mom_new << std::endl;

          const double mass = pdt_->particle(particle->pdgId()).ref().mass().value();
          const CLHEP::HepLorentzVector fourMom(mom_new, sqrt(mass*mass + mom_new.mag2()));

          output->push_back(GenParticle(particle->pdgId(),
                                        GenId::fromStepPointMCs,
                                        pos_new,
                                        fourMom,
                                        i->time(),
                                        i->properTime()
                                        ));

          history->addSingle(art::Ptr<GenParticle>(gpc_pid, output->size()-1, event.productGetter(gpc_pid)),
                             art::Ptr<StepPointMC>(makeMyHandle(ih), std::distance(ih->begin(), i))
                            );
        }
        else {
          mf::LogWarning warn("FromStepPointMCsRotateTarget");
          warn<<"WARNING: No particle data for pdgId = "<<particle->pdgId()<<"\n";
        }

      } // for (StepPointMCCollection entries)
    } // for(inputs_)

    event.put(std::move(output));
    event.put(std::move(history));
  }
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::FromStepPointMCsRotateTarget);
