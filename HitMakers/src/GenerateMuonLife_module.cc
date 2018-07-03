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

    explicit GenerateMuonLife(const fhicl::ParameterSet& pset);

    virtual void produce(art::Event& e) override;

  private:
    CLHEP::RandExponential  rexp_;
    double mean_;
    int  verbosityLevel_;
    bool readmap_;
    art::ProductToken<SimParticleTimeMap> const inmap_; // optional input map

    static double getMean(const  fhicl::ParameterSet& pset);
  };

  //================================================================
  GenerateMuonLife::GenerateMuonLife(const fhicl::ParameterSet& pset)
    : rexp_(createEngine( art::ServiceHandle<SeedService>()->getSeed() ))
    , mean_(getMean(pset))
    , verbosityLevel_(pset.get<int>("verbosityLevel", 0))
    , readmap_(pset.get<bool>("ReadMap",false))
    , inmap_{mayConsume<SimParticleTimeMap>(pset.get<art::InputTag>("InputTimeMap","null"))}
  {
    consumesMany<SimParticleCollection>();
    produces<SimParticleTimeMap>();
    if(verbosityLevel_ > 0) {
      std::cout<<"GenerateMuonLife initialized with meanLife = "<<mean_<<std::endl;
    }
  }

  //================================================================
  double GenerateMuonLife::getMean(const  fhicl::ParameterSet& pset) {
    double res(0);
    if(!pset.get_if_present("meanLife", res)) {
      // use the default
      GlobalConstantsHandle<PhysicsParams> phyPar;
      res = phyPar->getDecayTime();
    }
    return res;
  }

  //================================================================
  void GenerateMuonLife::produce(art::Event& event) {
  // optionally copy over input map
    SimParticleTimeMap* map(0);
    if(readmap_){
      auto inmap = event.getValidHandle(inmap_);
      map = new SimParticleTimeMap(*inmap.product());
    } else
      map = new SimParticleTimeMap;

    std::unique_ptr<SimParticleTimeMap> res(map);

    std::vector<art::Handle<SimParticleCollection> > colls;
    event.getManyByType(colls);

    // Generate and record offsets for all primaries
    for(const auto& ih : colls) {
      for(const auto& iter : *ih) {
        if(iter.second.isPrimary()) {
          art::Ptr<SimParticle> part(ih, iter.first.asUint());
	  // don't re-simulate if particle is already present.  This can happen if there is an input map
	  if(res->find(part) == res->end()){
	    if(part->genParticle()->generatorId() == GenId::StoppedParticleReactionGun
		|| part->genParticle()->generatorId() == GenId::conversionGun
		|| part->genParticle()->generatorId() == GenId::dioTail
	      )
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
