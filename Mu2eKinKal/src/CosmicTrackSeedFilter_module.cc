#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Principal/Handle.h"

#include "Offline/RecoDataProducts/inc/CosmicTrack.hh"
#include "Offline/RecoDataProducts/inc/CosmicTrackSeed.hh"
#include<TObject.h>
#include <iostream>

using namespace std;
namespace mu2e {

  class CosmicsTrackSeedFilter : public art::EDFilter {
    public:
      struct Config {
        using Name=fhicl::Name;
        using Comment=fhicl::Comment;
        fhicl::Atom<art::InputTag> cstag{Name("CosmicTrackSeedCollection"),Comment("tag for cosmic track seed collection")};
      };

      typedef art::EDFilter::Table<Config> Parameters;
      explicit CosmicsTrackSeedFilter(const Parameters& conf);
      virtual bool filter (art::Event& event) override;
      virtual ~CosmicsTrackSeedFilter() {}

    private:
      Config _conf;
      const CosmicTrackSeedCollection* _cscol;
      art::InputTag   _cstag;
      bool findData(const art::Event& evt);
  };

  CosmicsTrackSeedFilter::CosmicsTrackSeedFilter(const Parameters& conf):
    art::EDFilter{conf},
    _cstag(conf().cstag())
    {}


  bool CosmicsTrackSeedFilter::filter(art::Event& event) {
    bool pass = false;
    if(findData(event)) {
      pass = true;
    }
    return pass;
  }

  bool CosmicsTrackSeedFilter::findData(const art::Event& evt){
    _cscol = 0;
    auto csH = evt.getValidHandle<CosmicTrackSeedCollection>(_cstag);
    _cscol =csH.product();
    return  _cscol !=0;
  }
}
using mu2e::CosmicsTrackSeedFilter;
DEFINE_ART_MODULE(CosmicsTrackSeedFilter)
