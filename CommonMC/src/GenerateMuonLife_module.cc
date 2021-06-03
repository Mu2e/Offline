// Attach random, exponentially distributed time offsets to SimParticles
// made by the StoppedParticleReactionGun generator, and 0 offsets to
// particles from other origins.
//
// 20140211: conversionGun handling is added; this is a kludge to keep
// using the existing Kalman fit tuning and analysis code.
//
// Andrei Gaponenko, 2014

#include <string>
#include <memory>

#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "CLHEP/Random/RandExponential.h"

#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"

#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/SimParticleTimeMap.hh"
#include "MCDataProducts/inc/GenId.hh"
#include "SeedService/inc/SeedService.hh"
#include "Mu2eUtilities/inc/SimParticleCollectionPrinter.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/PhysicsParams.hh"

namespace mu2e {

  class GenerateMuonLife : public art::EDProducer {
  public:

    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::Atom<double> meanLife{
        Name("meanLife"),
          Comment("A negative value (default) means use the muon life time in the stopping target material\n"
                  "as given by the GlobalConstantsService. The particle life time can be overridden\n"
                  "by providing a positive number here, e.g. for simulating stopped pion daughters."
                  ),
          -1.
          };

      fhicl::Sequence<art::InputTag> InputTimeMaps {
        Name("InputTimeMaps"),
          Comment("Pre-exising time offsets that should be transferred to the output."),
          std::vector<art::InputTag>()
          };

      fhicl::Atom<int> verbosityLevel{ Name("verbosityLevel"), Comment("Levels 0, 1, and 11 increase the number of printouts.."), 0 };

      fhicl::Sequence<std::string> applyToGenIds {
        Name("applyToGenIds"),
          Comment("The whitelist mode: assign time offsets just to particles made by one of the\n"
                  "listed generators.\n")
          };
    };

    using Parameters = art::EDProducer::Table<Config>;
    explicit GenerateMuonLife(const Parameters& conf);

    virtual void produce(art::Event& e) override;

  private:
    CLHEP::RandExponential  rexp_;
    double mean_;
    int  verbosityLevel_;
    std::vector<art::ProductToken<SimParticleTimeMap> > inmaps_; // optional input maps

    typedef std::set<GenId::enum_type> GenIdSet;
    GenIdSet applyToGenIds_;
  };

  //================================================================
  GenerateMuonLife::GenerateMuonLife(const Parameters& conf)
    : EDProducer{conf}
    , rexp_(createEngine( art::ServiceHandle<SeedService>()->getSeed() ))
    , mean_(conf().meanLife())
    , verbosityLevel_(conf().verbosityLevel())
    {
      std::vector<art::InputTag> inmaps = conf().InputTimeMaps();
      for(auto const& tag : inmaps ){
        inmaps_.push_back(consumes<SimParticleTimeMap>(tag));
      }
      consumesMany<SimParticleCollection>();
      produces<SimParticleTimeMap>();

      if(mean_ <= 0.) {
        GlobalConstantsHandle<PhysicsParams> phyPar;
        mean_ = phyPar->getDecayTime();
      }
      if(verbosityLevel_ > 0) {
        std::cout<<"GenerateMuonLife initialized with meanLife = "<<mean_<<std::endl;
      }

      for(const auto i: conf().applyToGenIds()) {
        applyToGenIds_.insert(GenId::findByName(i).id());
      }
    }

  //================================================================
  void GenerateMuonLife::produce(art::Event& event) {
    std::unique_ptr<SimParticleTimeMap> res(new SimParticleTimeMap);
    // copy over input maps (if any)
    for(auto const& token : inmaps_) {
      auto inmap = event.getValidHandle(token);
      res->insert(inmap->begin(),inmap->end());
    }

    std::vector<art::Handle<SimParticleCollection> > colls = event.getMany<SimParticleCollection>();

    // Generate and record offsets for all primaries
    for(const auto& ih : colls) {
      for(const auto& iter : *ih) {
        if(iter.second.isPrimary()) {
          art::Ptr<SimParticle> part(ih, iter.first.asUint());
          // don't re-simulate if particle is already present.  This can happen if there is an input map
          if(res->find(part) == res->end()){
            const auto genId = part->genParticle()->generatorId();

            // do just explicitly listed GenIds
            bool apply= applyToGenIds_.find(genId.id()) != applyToGenIds_.end();
            if (!apply && genId.isConversion()) { // also want to anything that is isConversion
              apply = true;
            }
            if (apply)
              {
                (*res)[part] = rexp_.fire(mean_);
              }
            else
              {
                (*res)[part] = 0;
              }
          }
        }
      }
    }

    if(verbosityLevel_ > 10) {
      std::cout<<"GenerateMuonLife dump begin"<<std::endl;
      for(const auto& i : *res) {
        SimParticleCollectionPrinter::print(std::cout, *i.first);
        std::cout<<" => "<<i.second<<std::endl;
      }
      std::cout<<"GenerateMuonLife dump end"<<std::endl;
    }

    event.put(std::move(res));
  } // end GenerateMuonLife::produce.

  //================================================================
} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::GenerateMuonLife)
