// Attach random, flat distributed time offsets to SimParticles
// made by the Cosmic generator, and 0 offsets to
// particles from other origins.
//
// Yuri Oksuzian, 2019

#include <string>
#include <memory>

#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "CLHEP/Random/RandFlat.h"

#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"

#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/SimParticleTimeMap.hh"
#include "MCDataProducts/inc/GenId.hh"
#include "SeedService/inc/SeedService.hh"
#include "Mu2eUtilities/inc/SimParticleCollectionPrinter.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/PhysicsParams.hh"

namespace mu2e {

  class GenerateCosmicTimes : public art::EDProducer {

  public:

    explicit GenerateCosmicTimes(const fhicl::ParameterSet& pset);

    virtual void produce(art::Event& e) override;

  private:
    CLHEP::RandFlat _randflat;
    int verbosityLevel_, offsetToTracker_;
    double tmin_, tmax_;
    art::InputTag hitsInputTag_;
    std::vector<art::ProductToken<SimParticleTimeMap> > inmaps_; // optional input maps
  };

  //================================================================
  GenerateCosmicTimes::GenerateCosmicTimes(const fhicl::ParameterSet& pset)
    : art::EDProducer{pset}
    , _randflat(createEngine( art::ServiceHandle<SeedService>()->getSeed() ))
    , verbosityLevel_(pset.get<int>("verbosityLevel", 0))
    , offsetToTracker_(pset.get<bool>("offsetToTracker", true))
    , tmin_(pset.get<double>("tmin"))
    , tmax_(pset.get<double>("tmax"))
    , hitsInputTag_(pset.get<std::string>("hitsInputTag"))
  {
    std::vector<art::InputTag> inmaps = pset.get<std::vector<art::InputTag> >("InputTimeMaps",std::vector<art::InputTag>());
    for(auto const& tag : inmaps ){
      inmaps_.push_back(consumes<SimParticleTimeMap>(tag));
    }
    consumesMany<SimParticleCollection>();
    produces<SimParticleTimeMap>();
    if(verbosityLevel_ > 0) {
      std::cout<<"GenerateCosmicTimes initialized with range = ["<< tmin_ << ";"<< tmax_ << "]"<< std::endl;
    }
  }

  //================================================================
  void GenerateCosmicTimes::produce(art::Event& event) {
    std::unique_ptr<SimParticleTimeMap> res(new SimParticleTimeMap);
    // copy over input maps (if any)
    for(auto const& token : inmaps_) {
      auto inmap = event.getValidHandle(token);
      res->insert(inmap->begin(),inmap->end());
    }

    std::vector<art::Handle<SimParticleCollection> > colls;
    event.getManyByType(colls);

    art::Handle<std::vector<mu2e::StepPointMC>> spHndl;
    bool gotIt = event.getByLabel(hitsInputTag_, spHndl);
    
    double firstTrackerHit = 0;
    if(gotIt && offsetToTracker_ && spHndl->size() > 0){
      firstTrackerHit = FLT_MAX;
      std::vector<mu2e::StepPointMC> stepPoints = *spHndl;
      for (auto const& spmc : stepPoints) 
	if(spmc.time() < firstTrackerHit) firstTrackerHit = spmc.time();
    }

    // Generate and record offsets for all primaries
    for(const auto& ih : colls) {
      for(const auto& iter : *ih) {
        if(iter.second.isPrimary()) {
          art::Ptr<SimParticle> part(ih, iter.first.asUint());
	  // don't re-simulate if particle is already present.  This can happen if there is an input map
	  if(res->find(part) == res->end()){
	    if(part->genParticle()->generatorId() == GenId::cosmicCRY   ||
	       part->genParticle()->generatorId() == GenId::cosmicDYB   ||
	       part->genParticle()->generatorId() == GenId::cosmic )
	      {
		(*res)[part] = _randflat.fire(tmin_ - firstTrackerHit, tmax_ - firstTrackerHit);
		if(verbosityLevel_ > 0)
		  std::cout << tmin_ << " " << tmax_ << " " << firstTrackerHit << std::endl;
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
      std::cout<<"GenerateCosmicTimes dump begin"<<std::endl;
      for(const auto& i : *res) {
        SimParticleCollectionPrinter::print(std::cout, *i.first);
        std::cout<<" => "<<i.second<<std::endl;
      }
      std::cout<<"GenerateCosmicTimes dump end"<<std::endl;
    }

    event.put(std::move(res));
  } // end GenerateCosmicTimes::produce.

  //================================================================
} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::GenerateCosmicTimes)
