// Purpose: Event filter for RPC simulations
// author: S Middleton  2024
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"

// Mu2e includes.
#include "Offline/MCDataProducts/inc/StageParticle.hh"
#include "Offline/SeedService/inc/SeedService.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include <iostream>
#include <string>

#include "TTree.h"
using namespace std;
namespace mu2e {

  class PionFilter : public art::EDFilter {
    public:
      struct Config {
        using Name=fhicl::Name;
        using Comment=fhicl::Comment;
        fhicl::Atom<art::InputTag> SimToken{Name("SimParticleCollection"),Comment("")};
        fhicl::Atom<art::InputTag>SimTag{Name("SimParticleCollection"),Comment("SimTag")};
        fhicl::Atom<double> tmin{Name("tmin"),0};
        fhicl::Atom<double> tmax{Name("tmax"),1e6};
        fhicl::Atom<bool> isNull{Name("isNull"),true};
      };
      explicit PionFilter(const art::EDFilter::Table<Config>& config);
      virtual bool filter(art::Event& event) override;

    private:
      art::InputTag _SimToken;
      const SimParticleCollection* _SimCol;
      double tmin_;
      double tmax_;
      bool isNull_;
      TTree* genTree;
      Float_t _endglobaltime;
      Float_t _startglobaltime;
      Float_t _endpropertime;
      Float_t _startpropertime;
      double totalweight = 0;
      double selectedweight = 0;
  };

  PionFilter::PionFilter(const art::EDFilter::Table<Config>& config) :
     EDFilter{config}
    , _SimToken(config().SimToken())
    , tmin_{config().tmin()}
    , tmax_{config().tmax()}
    , isNull_{config().isNull()}
  {
      art::ServiceHandle<art::TFileService> tfs;
      genTree  = tfs->make<TTree>("GenAna", "GenAna");
      genTree->Branch("endglobaltime", &_endglobaltime, "endglobaltime/F");
      genTree->Branch("startglobaltime", &_startglobaltime, "startglobaltime/F");
      genTree->Branch("endpropertime", &_endpropertime, "endpropertime/F");
      genTree->Branch("startpropertime", &_startpropertime, "startpropertime/F");
  }

  bool PionFilter::filter(art::Event& evt) {
      if(isNull_) return true;
      bool passed = false;
      std::vector<art::Handle<SimParticleCollection>> vah = evt.getMany<SimParticleCollection>();
      for (auto const& ah : vah) { //always one collection
        for(const auto& aParticle : *ah){
          art::Ptr<SimParticle> pp(ah, aParticle.first.asUint());
          _endglobaltime = pp->endGlobalTime();
          _startglobaltime = pp->startGlobalTime();
          _endpropertime = pp->endProperTime();
          _startpropertime = pp->startProperTime();
          if( pp->stoppingCode() == ProcessCode::mu2eKillerVolume and abs(pp->pdgId())  == 211){
            const PhysicsParams& gc = *GlobalConstantsHandle<PhysicsParams>();
            totalweight += exp(-1*pp->endProperTime() / gc.getParticleLifetime(pp->pdgId()));
          }
          if( pp->stoppingCode() == ProcessCode::mu2eKillerVolume and (abs(pp->pdgId())  == 211 and _endglobaltime > tmin_ and _endglobaltime < tmax_ )){
            passed = true;
            const PhysicsParams& gc = *GlobalConstantsHandle<PhysicsParams>();
            selectedweight += exp(-1*pp->endProperTime() / gc.getParticleLifetime(pp->pdgId()));
            genTree->Fill();
          }
        }
      }
    //std::cout<<"Total weight for all stops "<<totalweight<<std::endl;
    //std::cout<<"Selected weight for chosen stops "<<selectedweight<<std::endl;
    return passed;
  }
}

using mu2e::PionFilter;
DEFINE_ART_MODULE(PionFilter)
