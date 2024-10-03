// Author: S Middleton
// Date: Dec 2019
// Purpose: For Developing Straight Track Trigger

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "Offline/RecoDataProducts/inc/TrkFitFlag.hh"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Sequence.h"

#include "Offline/RecoDataProducts/inc/CosmicTrackSeed.hh"
#include "Offline/RecoDataProducts/inc/TriggerInfo.hh"

#include <string>
#include <vector>
#include <iostream>
#include <memory>

using namespace CLHEP;

namespace mu2e
{
  class CosmicSeedFilter : public art::EDFilter
  {
  public:
    struct Config{
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<art::InputTag>   cosmictag {Name("CosmicTrackSeedCollection"),Comment("track collection")};
      fhicl::Sequence<std::string> goodcosmic{Name("cosmicseedFitFlag"),Comment("Required flags"),std::vector<std::string>{"HelixOK","HelixConverged"}};
      fhicl::Atom<unsigned>        minNStrawHits        {    Name("minNStrawHits"),                   Comment("minNStrawHits")};
      fhicl::Atom<float>           minFracStrawHits{Name("minFracStrawHits"), Comment("Minimum fraction of hits in time cluster on track")};
      fhicl::Atom<float>           maxTheta{Name("maxTheta"), Comment("Max track angle from vertical")};
      fhicl::Atom<unsigned>        minNPanelStereo{Name("minNPanelStereo"), Comment("Min number of 12 possible panel orientations on track")};
      fhicl::Atom<int>             debug     {Name("debugLevel"), Comment("set to 1 for debug prints")};
    };
    typedef art::EDFilter::Table<Config> Parameters;
    explicit CosmicSeedFilter(const Parameters& conf);

    virtual bool filter(art::Event& event) override;
    virtual bool endRun( art::Run& run )   override;

  private:

    art::InputTag   _cosmicTag;
    TrkFitFlag      _goodcosmic;
      unsigned      _minnhits;
      float         _minFracStrawHits;
      float         _maxTheta;
      unsigned      _minNPanelStereo;
    int             _debug;
    unsigned        _nevt, _npass;
  };

  CosmicSeedFilter::CosmicSeedFilter(const Parameters& conf) :
    art::EDFilter(conf),
    _cosmicTag (conf().cosmictag()),
    _goodcosmic (conf().goodcosmic()),
    _minnhits(conf().minNStrawHits()),
    _minFracStrawHits(conf().minFracStrawHits()),
    _maxTheta(conf().maxTheta()),
    _minNPanelStereo(conf().minNPanelStereo()),
    _debug (conf().debug()),
    _nevt(0), _npass(0)
  {
    produces<TriggerInfo>();
  }

  bool CosmicSeedFilter::filter(art::Event& evt){
    std::unique_ptr<TriggerInfo> triginfo(new TriggerInfo);
    ++_nevt;
    bool   retval(false);
    // find the collection
    auto cosH = evt.getValidHandle<CosmicTrackSeedCollection>(_cosmicTag);
    const CosmicTrackSeedCollection* coscol = cosH.product();
    for(auto icos = coscol->begin(); icos != coscol->end(); ++icos) {
      auto const& cs = *icos;

      if( cs.status().hasAllProperties(_goodcosmic) && cs.hits().size() > _minnhits &&
          (float)cs.hits().size() >= cs.timeCluster()->nStrawHits()*_minFracStrawHits &&
          acos(fabs(cs.track().FitEquation.Dir.Unit().y())) <= _maxTheta){
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
        if (panelcount >= _minNPanelStereo)
        {


          ++_npass;
          retval = true;
          size_t index = std::distance(coscol->begin(),icos);
          triginfo->_cosmics.push_back(art::Ptr<CosmicTrackSeed>(cosH,index));
          if(_debug > 1){
            std::cout << moduleDescription().moduleLabel() << " passed event " << evt.id() << std::endl;
          }
          break;
        }
      }
    }
    evt.put(std::move(triginfo));
    return retval;
  }

  bool CosmicSeedFilter::endRun( art::Run& run ) {
    if(_debug > 0 && _nevt > 0){
      std::cout << moduleDescription().moduleLabel() << " passed " <<  _npass << " events out of " << _nevt << " for a ratio of " << float(_npass)/float(_nevt) << std::endl;
    }
    return true;
  }
}
using mu2e::CosmicSeedFilter;
DEFINE_ART_MODULE(CosmicSeedFilter)
