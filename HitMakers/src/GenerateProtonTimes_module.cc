// Draw time values from the proton time profile distribution and
// associate SimParticles with time offsets of their primary protons.
//
// Andrei Gaponenko, 2013

#include <string>
#include <memory>

#include "fhiclcpp/ParameterSet.h"
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
#include "SeedService/inc/SeedService.hh"
#include "Mu2eUtilities/inc/ProtonPulseRandPDF.hh"
#include "Mu2eUtilities/inc/SimParticleCollectionPrinter.hh"

namespace mu2e {

  class GenerateProtonTimes : public art::EDProducer {

  public:

    explicit GenerateProtonTimes(const fhicl::ParameterSet& pset);

    virtual void beginRun(art::Run&   r) override;
    virtual void produce (art::Event& e) override;

  private:
    art::RandomNumberGenerator::base_engine_t& engine_;
    fhicl::ParameterSet protonPset_;
    int  verbosityLevel_;

    typedef std::set<GenId::enum_type> GenIdSet;
    GenIdSet ignoredGenIds_;
    GenIdSet applyToGenIds_;

    std::unique_ptr<ProtonPulseRandPDF>  protonPulse_;

    std::string listStream( const GenIdSet& vsList );

  };

  //================================================================
  GenerateProtonTimes::GenerateProtonTimes(fhicl::ParameterSet const& pset)
    : engine_(createEngine(art::ServiceHandle<SeedService>()->getSeed()) )
    , protonPset_( pset.get<fhicl::ParameterSet>("randPDFparameters", fhicl::ParameterSet() ) )
    , verbosityLevel_(pset.get<int>("verbosityLevel", 0)) 
  {
    produces<SimParticleTimeMap>();

    typedef std::vector<std::string> VS;

    const auto ig(pset.get<VS>("ignoredGenIds", VS{"cosmicToy", "cosmicDYB", "cosmic"}));
    for(const auto i: ig) {
      ignoredGenIds_.insert(GenId::findByName(i).id());
    }

    const auto ia(pset.get<VS>("applyToGenIds", VS()));
    for(const auto i: ia) {
      applyToGenIds_.insert(GenId::findByName(i).id());
    }

    if(!applyToGenIds_.empty() && !ignoredGenIds_.empty()) {
      throw cet::exception("BADCONFIG")
        <<"If applyToGenIds is set, ignoredGenIds should be empty.  Got: applyToGenIds = [ "
        <<listStream( applyToGenIds_ )<<" ], ignoredGenIds = [ "<<listStream( ignoredGenIds_ )<<"]\n";
    }

  }

  //================================================================
  void GenerateProtonTimes::beginRun(art::Run& run) {
    protonPulse_.reset( new ProtonPulseRandPDF( engine_, protonPset_ ) );

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

    std::vector<art::Handle<SimParticleCollection> > colls;
    event.getManyByType(colls);

    // Generate and record offsets for all primaries
    for(const auto& ih : colls) {
      for(const auto& iter : *ih) {
        if(iter.second.isPrimary()) {
          art::Ptr<SimParticle> part(ih, iter.first.asUint());
          const auto genId = part->genParticle()->generatorId();

          const bool apply= applyToGenIds_.empty() ?
            // do all particles, except the explicitly vetoed ones
            (ignoredGenIds_.find(genId.id()) == ignoredGenIds_.end())
            // do just explicitly listed GenIds
            : (applyToGenIds_.find(genId.id()) != applyToGenIds_.end());

          (*res)[part] = apply ? protonPulse_->fire() : 0.;
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
