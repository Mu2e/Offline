// Purpose: Event filter for RPC simulations
// author: S Middleton  2024
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"

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
#include "TH1F.h"
#include "TTree.h"
using namespace std;
namespace mu2e {

  class PionFilter : public art::EDFilter {
    public:
      struct Config {
        using Name=fhicl::Name;
        using Comment=fhicl::Comment;
        fhicl::Atom<int> diagLevel{Name("diagLevel")};
        fhicl::Atom<double> tmin{Name("tmin")};
        fhicl::Atom<double> tmax{Name("tmax")};
        fhicl::Atom<int> processCode{Name("processCode")};
        fhicl::Atom<bool> isNull{Name("isNull")};
      };
      explicit PionFilter(const art::EDFilter::Table<Config>& config);
      virtual bool filter(art::Event& event) override;
      virtual void beginJob() override;
      virtual void endJob() override;

    private:
      const SimParticleCollection* SimCol_;
      int diagLevel_;
      double tmin_;
      double tmax_;
      int processCode_;
      bool isNull_;
      float _totalweight = 0;
      float _selectedweight = 0;
      int _ntot = 0;
      int _nsel = 0;
  };

  PionFilter::PionFilter(const art::EDFilter::Table<Config>& config) :
     EDFilter{config}
    , diagLevel_{config().diagLevel()}
    , tmin_{config().tmin()}
    , tmax_{config().tmax()}
    , processCode_{config().processCode()}
    , isNull_{config().isNull()}
  {
  }

  void PionFilter::beginJob(){
      art::ServiceHandle<art::TFileService> tfs;

}

  bool PionFilter::filter(art::Event& evt) {
      if(isNull_) return true;
      bool passed = false;
      std::vector<art::Handle<SimParticleCollection>> vah = evt.getMany<SimParticleCollection>();
      for (auto const& ah : vah) { //always one collection
        for(const auto& aParticle : *ah){
          art::Ptr<SimParticle> pp(ah, aParticle.first.asUint());

          float _endglobaltime = pp->endGlobalTime();
          if( pp->stoppingCode() == processCode_ and std::abs(pp->pdgId()) == PDGCode::pi_plus){
            const PhysicsParams& gc = *GlobalConstantsHandle<PhysicsParams>();
            _totalweight += exp(-1*pp->endProperTime() / gc.getParticleLifetime(pp->pdgId()));
            _ntot += 1;
          }
          if( pp->stoppingCode() == processCode_ and std::abs(pp->pdgId()) == PDGCode::pi_plus and _endglobaltime > tmin_ and _endglobaltime < tmax_ ){
            passed = true;
            const PhysicsParams& gc = *GlobalConstantsHandle<PhysicsParams>();
            _selectedweight += exp(-1*pp->endProperTime() / gc.getParticleLifetime(pp->pdgId()));
            _nsel += 1;
          }
        }
      }
    return passed;
  }

  void PionFilter::endJob(){
     if(diagLevel_ > 0 ){
      std::cout<<"Total weight for all stops "<<_totalweight<<std::endl;
      std::cout<<"Total stops "<<_ntot<<std::endl;
      std::cout<<"Selected weight for chosen stops "<<_selectedweight<<std::endl;
      std::cout<<"Selected stops "<<_nsel<<std::endl;
    }
  }
}

using mu2e::PionFilter;
DEFINE_ART_MODULE(PionFilter)
