#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
// mu2e
// data
#include "Offline/RecoDataProducts/inc/CosmicTrackSeed.hh"
#include "Offline/RecoDataProducts/inc/TriggerInfo.hh"
// c++
#include <iostream>
#include <memory>


namespace mu2e
{
  class CosmicTrackSeedFilter : public art::EDFilter
  {
    public:
      struct Config{
        using Name    = fhicl::Name;
        using Comment = fhicl::Comment;
        fhicl::Atom<art::InputTag>      csCollection{    Name("CosmicTrackSeedCollection"),      Comment("CosmicTrackSeedCollection label") };
        fhicl::Atom<unsigned>           minNStrawHits        {    Name("minNStrawHits"),                   Comment("minNStrawHits")};
        fhicl::Atom<float>              minFracStrawHits{Name("minFracStrawHits"), Comment("Minimum fraction of hits in time cluster on track")};
        fhicl::Atom<float>              maxTheta{Name("maxTheta"), Comment("Max track angle from vertical")};
        fhicl::Atom<unsigned>           minNPanelStereo{Name("minNPanelStereo"), Comment("Min number of 12 possible panel orientations on track")};
        fhicl::Atom<int>                debugLevel           {    Name("debugLevel"),                 Comment("Debug"),0 };
        fhicl::Atom<int>                noFilter             {    Name("noFilter"),                 Comment("Don't filter anything"),0 };
      };

      using Parameters = art::EDFilter::Table<Config>;

      explicit CosmicTrackSeedFilter(const Parameters& config);

    private:
      bool filter(art::Event& event) override;
      bool endRun(art::Run& run ) override;

      art::InputTag _csTag;
      unsigned      _minnhits;
      float         _minFracStrawHits;
      float         _maxTheta;
      unsigned      _minNPanelStereo;
      int           _debug;
      // counters
      unsigned      _nevt, _npass;
      int           _noFilter;
  };

  CosmicTrackSeedFilter::CosmicTrackSeedFilter(const Parameters& conf)
    : art::EDFilter{conf},
    _csTag   (conf().csCollection()),
    _minnhits(conf().minNStrawHits()),
    _minFracStrawHits(conf().minFracStrawHits()),
    _maxTheta(conf().maxTheta()),
    _minNPanelStereo(conf().minNPanelStereo()),
    _debug   (conf().debugLevel()),
    _nevt    (0),
    _npass   (0),
    _noFilter(conf().noFilter())
    {
      produces<TriggerInfo>();
    }

  bool CosmicTrackSeedFilter::filter(art::Event& evt){
    // create output
    std::unique_ptr<TriggerInfo> triginfo(new TriggerInfo);
    ++_nevt;
    bool retval(false); // preset to fail
    // find the collection
    auto csH = evt.getValidHandle<CosmicTrackSeedCollection>(_csTag);
    const CosmicTrackSeedCollection* cscol = csH.product();
    // loop over the collection: if any pass the selection, pass this event
    for(auto ics = cscol->begin();ics != cscol->end(); ++ics) {
      auto const& cs = *ics;
      if(_debug > 2){
        std::cout << moduleDescription().moduleLabel() << " nStrawHits = " << cs.hits().size() << " t0 = " << cs.t0().t0() << std::endl;
      }

      std::vector<bool> panel(12,false);
      for (size_t k=0;k<cs.hits().size();k++){
        auto sid = cs.hits()[k].strawId();
        if ((sid.plane()%4) == 0 || (sid.plane()%4) == 3){
          panel[sid.panel()] = true;
        }else{
          panel[sid.panel()+6] = true;
        }
      }

      unsigned panelcount = 0;
      for (size_t i=0;i<12;i++){
        if (panel[i])
          panelcount++;
      }
      if (cs.hits().size() >= _minnhits &&
          (float) cs.hits().size()  >= cs.timeCluster()->nStrawHits()*_minFracStrawHits &&
          acos(fabs(cs.track().FitEquation.Dir.Unit().y())) <= _maxTheta &&
          panelcount >= _minNPanelStereo)
      {
        retval = true;
        ++_npass;
        // Fill the trigger info object
        // associate to the hit cluster which triggers.  Note there may be other hit clusters which also pass the filter
        // but filtering is by event!
        size_t index = std::distance(cscol->begin(),ics);
        triginfo->_cosmics.push_back(art::Ptr<CosmicTrackSeed>(csH,index));

        if(_debug > 1){
          std::cout << moduleDescription().moduleLabel() << " passed event " << evt.id() << std::endl;
        }
      }
    }
    evt.put(std::move(triginfo));
    if (_noFilter != 1){
      return retval;
    }else {
      return true;
    }
  }

  bool CosmicTrackSeedFilter::endRun( art::Run& run ) {
    if(_debug > 0 && _nevt > 0){
      std::cout << moduleDescription().moduleLabel() << " passed " << _npass << " events out of " << _nevt << " for a ratio of " << float(_npass)/float(_nevt) << std::endl;
    }
    return true;
  }
}
using mu2e::CosmicTrackSeedFilter;
DEFINE_ART_MODULE(CosmicTrackSeedFilter)
