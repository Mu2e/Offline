// Attach random, flat distributed time offsets to SimParticles
// made by the Cosmic generator, and 0 offsets to
// particles from other origins.  If detector steps are provided,
// use the earliest to set the minimum time
//
// Yuri Oksuzian, 2019

#include <string>
#include <memory>
#include <cfloat>

#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "CLHEP/Random/RandFlat.h"

#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"

#include "DataProducts/inc/EventWindowMarker.hh"
#include "MCDataProducts/inc/StrawGasStep.hh"
#include "MCDataProducts/inc/CaloShowerStep.hh"
#include "MCDataProducts/inc/CrvStep.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/SimParticleTimeMap.hh"
#include "MCDataProducts/inc/GenId.hh"
#include "SeedService/inc/SeedService.hh"
#include "Mu2eUtilities/inc/SimParticleCollectionPrinter.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "ProditionsService/inc/ProditionsHandle.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "TrackerConditions/inc/StrawElectronics.hh"

namespace mu2e {

  class GenerateCosmicTimes : public art::EDProducer {

  public:
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::Atom<int> verbosityLevel{ Name("verbosityLevel"), Comment("Verbosity Level"),0 };
      fhicl::Atom<double> tBuff{ Name("TimeBuffer"), Comment("Time Offset buffer, to account for propagation and digitization delays (ns)"),150.0 };
      fhicl::Atom<double> genTmin{ Name("TimeMin"), Comment("TEST: Minimum generation time  (ns)"), -1};
      fhicl::Atom<double> genTmax{ Name("TimeMax"), Comment("TEST: Maximum generation time  (ns)"), -1};
      fhicl::Sequence<art::InputTag> trkSteps { Name("StrawGasSteps"), Comment("StrawGasStep collections") };
      fhicl::Sequence<art::InputTag> caloSteps { Name("CaloShowerSteps"), Comment("CaloShowerStep collections") };
      fhicl::Sequence<art::InputTag> crvSteps { Name("CrvSteps"), Comment("CrvStep collections") };
      fhicl::Sequence<art::InputTag> inMaps { Name("InputTimeMaps"), Comment("Input time maps to append to"), std::vector<art::InputTag>{} };
      fhicl::Atom<art::InputTag> ewMarkerTag{ Name("EventWindowMarker"), Comment("EventWindowMarker producer"),"EWMProducer" };
    };

    using Parameters = art::EDProducer::Table<Config>;

    explicit GenerateCosmicTimes(const Parameters& conf);

    virtual void produce(art::Event& e) override;

  private:
    CLHEP::RandFlat _randflat;
    int verbosityLevel_;
    double tbuff_; // time buffer to use when defining the event window
    double genTmin_; // TEST: Minimum generation time
    double genTmax_; // TEST: Maximum generation time
    std::vector<art::InputTag> trkStepCols_, caloStepCols_, crvStepCols_;
    art::InputTag ewMarkerTag_; 
    std::vector<art::ProductToken<SimParticleTimeMap> > inmaps_; // optional input maps
    ProditionsHandle<StrawElectronics> strawele_h_;
  };

  //================================================================
  GenerateCosmicTimes::GenerateCosmicTimes(const Parameters& conf) : art::EDProducer{conf}
  , _randflat(createEngine( art::ServiceHandle<SeedService>()->getSeed() ))
    , verbosityLevel_(conf().verbosityLevel())
    , tbuff_(conf().tBuff())
    , genTmin_(conf().genTmin())
    , genTmax_(conf().genTmax())
    , ewMarkerTag_(conf().ewMarkerTag())
  {
    for(const auto& trktag : conf().trkSteps()) { trkStepCols_.emplace_back(trktag); consumes<StrawGasStepCollection>(trktag); }
    for(const auto& calotag : conf().caloSteps()) { caloStepCols_.emplace_back(calotag); consumes<CaloShowerStepCollection>(calotag); }
    for(const auto& crvtag : conf().crvSteps()) { crvStepCols_.emplace_back(crvtag);  consumes<CrvStepCollection>(crvtag); }

    for(auto const& tag : conf().inMaps() ){ inmaps_.push_back(consumes<SimParticleTimeMap>(tag)); consumes<SimParticleTimeMap>(tag); }
    consumesMany<SimParticleCollection>();
    consumes<EventWindowMarker>(ewMarkerTag_);
    produces<SimParticleTimeMap>();

  }

  //================================================================
  void GenerateCosmicTimes::produce(art::Event& event) {
  // create output
    std::unique_ptr<SimParticleTimeMap> res(new SimParticleTimeMap);
    // copy over input maps (if any)
    for(auto const& token : inmaps_) {
      auto inmap = event.getValidHandle(token);
      res->insert(inmap->begin(),inmap->end());
    }
    // find EventWindowMarker: this defines the event length
    art::Handle<EventWindowMarker> ewMarkerHandle;
    event.getByLabel(ewMarkerTag_, ewMarkerHandle);
    const EventWindowMarker& ewMarker(*ewMarkerHandle);
      
    ConditionsHandle<AcceleratorParams> accPar("ignored");
    auto mbtime = accPar->deBuncherPeriod;

// find the earliest step.
// Use this to define the earliest time, to improves the generation efficiency
    double tearly(FLT_MAX);
    for(const auto& trkcoltag : trkStepCols_) {
      auto sgscolH = event.getValidHandle<StrawGasStepCollection>(trkcoltag);
      for(const auto& sgs : *sgscolH ) {
	tearly = std::min(tearly,sgs.time());
      }
    }
 
    for(const auto& calocoltag : caloStepCols_) {
      auto csscolH = event.getValidHandle<CaloShowerStepCollection>(calocoltag);
      for(const auto& css : *csscolH ) {
	tearly = std::min(tearly,css.time());
      }
    }

    for(const auto& crvcoltag : crvStepCols_) {
      auto crvscolH = event.getValidHandle<CrvStepCollection>(crvcoltag);
      for(const auto& crvs : *crvscolH ) {
	tearly = std::min(tearly,crvs.startTime());	
      }
    }
    if(tearly >= FLT_MAX)tearly = 0.0;

// use a buffer to set the earliest time offset
    double tmin = -tearly - tbuff_;

// adjust start time according to digitization window start.
// this assumes the delayed digitization start WILL be applied also for offspill in the digitizers FIXME!
    StrawElectronics const& strawele = strawele_h_.get(event.id());
    tmin += strawele.digitizationStart(); // flash end should be a Mu2e global quantity FIXME!
    double tmax = ewMarker.eventLength() - tearly;
    if(ewMarker.spillType() == EventWindowMarker::onspill){
      tmax = mbtime - tearly;
    }

    if(genTmin_ > 0)
      tmin = genTmin_;
    if(genTmax_ > 0)
      tmax = genTmax_;

    if(verbosityLevel_ > 1) {
      std::cout << "DetectorStep early time = " << tearly << std::endl;
      std::cout<<"GenerateCosmicTimes time range = ["<< tmin << ","<< tmax << "]"<< std::endl;
    }



    // Generate and record offsets for all primaries
    std::vector<art::Handle<SimParticleCollection> > colls = event.getMany<SimParticleCollection>();
    for(const auto& ih : colls) {
      for(const auto& iter : *ih) {
        if(iter.second.isPrimary()) {
          art::Ptr<SimParticle> part(ih, iter.first.asUint());
	  // don't re-simulate if particle is already present.  This can happen if there is an input map
	  if(res->find(part) == res->end()){
	    if(part->genParticle()->generatorId() == GenId::cosmicCRY   ||
	       part->genParticle()->generatorId() == GenId::cosmicDYB   ||
	       part->genParticle()->generatorId() == GenId::cosmic      ||
               part->genParticle()->generatorId() == GenId::cosmicCORSIKA)
	      {
		(*res)[part] = _randflat.fire(tmin, tmax);
		if(verbosityLevel_ > 1)
		  std::cout << "Cosmic particle " << part->genParticle()->generatorId() << " given time " << (*res)[part] << std::endl;
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
