// Draw time values from the proton time profile distribution and
// associate SimParticles with time offsets of their primary protons.
//
// Andrei Gaponenko, 2013

#include <string>
#include <memory>

#include "fhiclcpp/types/OptionalTable.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "CLHEP/Random/RandGaussQ.h"

#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "canvas/Utilities/InputTag.h"

#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/SimParticleTimeMap.hh"
#include "MCDataProducts/inc/FixedTimeMap.hh"
#include "SeedService/inc/SeedService.hh"
#include "Mu2eUtilities/inc/ProtonPulseRandPDF.hh"
#include "Mu2eUtilities/inc/SimParticleCollectionPrinter.hh"

namespace mu2e {

  class GenerateProtonTimes : public art::EDProducer {

  public:

    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::OptionalTable<ProtonPulseRandPDF::Config> randPDFparameters { Name("randPDFparameters") };

      fhicl::Sequence<std::string> ignoredGenIds {
        Name("ignoredGenIds"),
          Comment("A non-empty list means: generate time offsets for all particles\n"
                  "except those produced by the listed generators.   If this list is emtpy,\n"
                  "the applyToGenIds parameter below is enabled."
                  ) };

      fhicl::Sequence<std::string> applyToGenIds {
        Name("applyToGenIds"),
          Comment("The whitelist mode: assign time offsets just to particles made by one of the\n"
                  "listed generators. This setting is only active in the case ignoredGenIds is emtpy.\n"
                  ),
       [this](){ return ignoredGenIds().empty();}
      };

      fhicl::Sequence<art::InputTag> InputTimeMaps {
        Name("InputTimeMaps"),
          Comment("Pre-exising time offsets that should be transferred to the output."),
          std::vector<art::InputTag>()
          };

      fhicl::Atom<art::InputTag> FixedModule {
        Name("FixedModule"),
          Comment("Input tag of a FixedTimeMap from RPC generation, if any."),
          art::InputTag()
          };

      fhicl::Atom<int> verbosityLevel{ Name("verbosityLevel"), Comment("Levels 0, 1, 3, and 11 increase the number of printouts.."), 0 };
    };

    using Parameters = art::EDProducer::Table<Config>;
    explicit GenerateProtonTimes(const Parameters& conf);

    virtual void beginRun(art::Run&   r) override;
    virtual void produce (art::Event& e) override;

  private:
    art::RandomNumberGenerator::base_engine_t& engine_;
    ProtonPulseRandPDF::Config protonPulseConf_;
    int  verbosityLevel_;
    art::InputTag fixedTime_;

    typedef std::set<GenId::enum_type> GenIdSet;
    GenIdSet ignoredGenIds_;
    GenIdSet applyToGenIds_;

    std::unique_ptr<ProtonPulseRandPDF>  protonPulse_;

    std::string listStream( const GenIdSet& vsList );
    std::vector<art::ProductToken<SimParticleTimeMap> > inmaps_; // optional input maps

  };

  //================================================================
  GenerateProtonTimes::GenerateProtonTimes(const Parameters& conf)
    : EDProducer{conf}
    , engine_(createEngine(art::ServiceHandle<SeedService>()->getSeed()) )
    , verbosityLevel_(conf().verbosityLevel())
    , fixedTime_(conf().FixedModule())
  {
      // require either ignore or applyto be non-null

    if(conf().ignoredGenIds().empty() && conf().applyToGenIds().empty() )
      throw cet::exception("Simulation")<<"No inclusion or exclusion GenIds specified" << std::endl;

    std::vector<art::InputTag> inmaps = conf().InputTimeMaps();
    for(auto const& tag : inmaps ){
      inmaps_.push_back(consumes<SimParticleTimeMap>(tag));
    }
    consumesMany<SimParticleCollection>();
    consumesMany<FixedTimeMap>();
    produces<SimParticleTimeMap>();

    typedef std::vector<std::string> VS;

    for(const auto i: conf().ignoredGenIds()) {
      ignoredGenIds_.insert(GenId::findByName(i).id());
    }

    if(ignoredGenIds_.empty()) {
      for(const auto i: conf().applyToGenIds()) {
        applyToGenIds_.insert(GenId::findByName(i).id());
      }
    }

    conf().randPDFparameters(protonPulseConf_);
  }

  //================================================================
  void GenerateProtonTimes::beginRun(art::Run& run) {
    protonPulse_.reset( new ProtonPulseRandPDF( engine_, protonPulseConf_ ) );

    if(verbosityLevel_ > 0) {
      if(applyToGenIds_.empty()) {
        mf::LogInfo("Info")<<"pulseType = "<<protonPulse_->pulseType() <<", ignoring genIds [ "<< listStream( ignoredGenIds_ ) <<" ]\n";
      }
      else {
        mf::LogInfo("Info")<<"pulseType = "<<protonPulse_->pulseType() <<", applying to genIds [ "<< listStream( applyToGenIds_ ) <<" ]\n";
      }
    }

    if ( verbosityLevel_ > 10 ) {
      std::ostringstream timeSpectrum;
      std::cout << " Size of proton pulse: " << protonPulse_->getTimes().size() << std::endl;
      for ( std::size_t i(0) ; i < protonPulse_->getTimes().size(); i++ ) {
        timeSpectrum << "   POT time: "
                     << protonPulse_->getTimes().at(i)
                     << "     "
                     << protonPulse_->getSpectrum().at(i) << "\n";
      }
      mf::LogInfo("Info") << "Longitudinal POT time distribution\n" << timeSpectrum.str();
    }
  }

  //================================================================
  void GenerateProtonTimes::produce(art::Event& event) {
    std::unique_ptr<SimParticleTimeMap> res(new SimParticleTimeMap);
    // copy over input maps (if any)
    for(auto const& token : inmaps_) {
      auto inmap = event.getValidHandle(token);
      res->insert(inmap->begin(),inmap->end());
    }

    std::vector<art::Handle<SimParticleCollection> > colls = event.getMany<SimParticleCollection>();
    art::Handle<FixedTimeMap> ftmHandle;
    if (!fixedTime_.empty())
      event.getByLabel(fixedTime_, ftmHandle);

    // Generate and record offsets for all primaries
    for(const auto& ih : colls) {
      for(const auto& iter : *ih) {
        if(iter.second.isPrimary()) {
          art::Ptr<SimParticle> part(ih, iter.first.asUint());
          // don't re-simulate if particle is already present.  This can happen if there is an input map
          if(res->find(part) == res->end()){
            const auto genId = part->genParticle()->generatorId();

            const bool apply= applyToGenIds_.empty() ?
              // do all particles, except the explicitly vetoed ones
              (ignoredGenIds_.find(genId.id()) == ignoredGenIds_.end())
              // do just explicitly listed GenIds
              : (applyToGenIds_.find(genId.id()) != applyToGenIds_.end());
	    if(verbosityLevel_ > 2){
	      if(apply)
		std::cout << "Applying proton time to genId " << genId << std::endl;
	      else
		std::cout << "Proton time NOT applied to genId " << genId << std::endl;
	    }
	      
            (*res)[part] = apply ? protonPulse_->fire() : 0.;
            if (!fixedTime_.empty()){
              (*res)[part] = apply ? ftmHandle->time() : 0;
            }
          } else if(verbosityLevel_ > 2) {
            std::cout << "Found existing particle in map" << std::endl;
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
  std::string GenerateProtonTimes::listStream( const GenIdSet& vsList ) {
    std::ostringstream osList;
    for(const auto i: vsList ) {
      osList<<i<<", ";
    }
    return osList.str();
  }

  //================================================================
} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::GenerateProtonTimes)
