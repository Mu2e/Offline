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
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include <iostream>
#include <string>

using namespace std;
namespace mu2e {

  class PionFilter : public art::EDFilter {
    public:
      struct Config {
        using Name=fhicl::Name;
        using Comment=fhicl::Comment;
        fhicl::Atom<int> diagLevel{Name("diagLevel"),0};
        fhicl::Atom<double> tmin{Name("tmin"),0};
        fhicl::Atom<double> tmax{Name("tmax"),1e6};
        fhicl::Atom<bool> isNull{Name("isNull"),true};
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
    , isNull_{config().isNull()}
  {
  }
  void PionFilter::beginJob(){}

  bool PionFilter::filter(art::Event& evt) {

      if(isNull_) return true;
      bool passed = false;
      std::vector<art::Handle<SimParticleCollection>> vah = evt.getMany<SimParticleCollection>();
      for (auto const& ah : vah) { //always one collection
        for(const auto& aParticle : *ah){
          art::Ptr<SimParticle> pp(ah, aParticle.first.asUint());
          float _endglobaltime = pp->endGlobalTime();
          if( pp->stoppingCode() == ProcessCode::hBertiniCaptureAtRest and std::abs(pp->pdgId()) == PDGCode::pi_plus){
            const PhysicsParams& gc = *GlobalConstantsHandle<PhysicsParams>();
            _totalweight += exp(-1*pp->endProperTime() / gc.getParticleLifetime(pp->pdgId()));
            _ntot += 1;
          }
          if( pp->stoppingCode() == ProcessCode::hBertiniCaptureAtRest and (std::abs(pp->pdgId()) == PDGCode::pi_plus and _endglobaltime > tmin_ and _endglobaltime < tmax_ )){
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
