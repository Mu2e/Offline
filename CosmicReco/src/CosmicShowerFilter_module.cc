//
//  Cut events where pattern of hits in a panel suggests showering
//  or a non-physical event
//  original author: Richie Bonventre (LBNL)
//

// Mu2e includes.
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/TrackerConditions/inc/TrackerStatus.hh"
// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Persistency/Common/Ptr.h"
// Other includes
#include "messagefacility/MessageLogger/MessageLogger.h"
// C++ includes
#include <iostream>
#include <memory>
#include <algorithm>
#include <utility>

using namespace std;

namespace mu2e {
  class CosmicShowerFilter : public art::EDProducer {
    public:
      struct Config{
        using Name=fhicl::Name;
        using Comment=fhicl::Comment;
        fhicl::OptionalAtom<int> maxnpanel {Name("maxNPanel"), Comment("maximum number of panels for shower cut")};
        fhicl::Atom<int> maxcrossinggap {Name("maxCrossingGap"), Comment("maximum missing straws when crossing layers")};
        fhicl::Atom<int> maxsamegap {Name("maxSameGap"), Comment("maximum missing straws on same layers"),1};
        fhicl::Atom<int> maxtotalsamegap {Name("maxTotalSameGap"), Comment("maximum total missing straws on same layers"),1};
        fhicl::Atom<bool> cutsinglelayer {Name("cutSingleLayer"), Comment("remove events on only a single layer of a panel"),false};
        fhicl::Atom<art::InputTag> tcToken{Name("TimeClusterCollection"),Comment("tag for time cluster collection")};
        fhicl::Atom<art::InputTag> chToken{Name("ComboHitCollection"),Comment("tag for combo hit collection")};
      };
      typedef art::EDProducer::Table<Config> Parameters;
      explicit CosmicShowerFilter(const Parameters& conf);

      void produce(art::Event& e) override;

    private:
      bool _hasmaxnpanel;
      int _maxnpanel;
      int _maxCrossingGap;
      int _maxSameGap;
      int _maxTotalSameGap;
      bool _cutsinglelayer;
      art::InputTag _tcToken;
      art::InputTag _chToken;

      ProditionsHandle<TrackerStatus> _trackerStatus_h;

  };

  CosmicShowerFilter::CosmicShowerFilter(const Parameters& conf) :
    art::EDProducer(conf),
    _hasmaxnpanel (false),
    _maxnpanel (0),
    _maxCrossingGap (conf().maxcrossinggap()),
    _maxSameGap (conf().maxsamegap()),
    _maxTotalSameGap (conf().maxtotalsamegap()),
    _cutsinglelayer (conf().cutsinglelayer()),
    _tcToken (conf().tcToken()),
    _chToken (conf().chToken())
  {
    consumes<TimeClusterCollection>(_tcToken);
    consumes<ComboHitCollection>(_chToken);
    produces<TimeClusterCollection>();
    _hasmaxnpanel = conf().maxnpanel(_maxnpanel);
  }


  //--------------------------------------------------------------------------------------------------------------
  void CosmicShowerFilter::produce(art::Event & event ){

    TrackerStatus const& trackerStatus = _trackerStatus_h.get(event.id());

    std::unique_ptr<TimeClusterCollection> tccol_filter(new TimeClusterCollection());

    auto const& tcH = event.getValidHandle<TimeClusterCollection>(_tcToken);
    const TimeClusterCollection& tccol(*tcH);
    auto const& chH = event.getValidHandle<ComboHitCollection>(_chToken);
    const ComboHitCollection& chcol(*chH);

    for (size_t index=0;index< tccol.size();++index) {
      const auto& tclust = tccol[index];

      // get the individual straw level hits
      ComboHitCollection shcol;
      std::vector<ComboHitCollection::const_iterator> chids;
      chcol.fillComboHits(event, tclust._strawHitIdxs, chids);
      for (auto const& it : chids){
        shcol.push_back(it[0]);
      }

      std::vector<std::vector<uint16_t> > straws(216,std::vector<uint16_t>());
      std::vector<uint16_t> panels;

      bool singlelayer = true;
      int singlelayer_panel = shcol[0].strawId().uniquePanel();
      int singlelayer_layer = shcol[0].strawId().layer();

      for(size_t ich = 0;ich < shcol.size(); ++ich){
        ComboHit const& chit = shcol[ich];
        if (chit.strawId().uniquePanel() != singlelayer_panel || chit.strawId().layer() != singlelayer_layer){
          singlelayer = false;
        }
        uint16_t panelid = chit.strawId().uniquePanel();
        if (std::find(panels.begin(),panels.end(),panelid) == panels.end())
          panels.push_back(panelid);
        straws[panelid].push_back(chit.strawId().getStraw());
      }
      // filter if too many panels in this cluster
      if (_hasmaxnpanel && (int) panels.size() > _maxnpanel)
        continue;
      if (singlelayer && _cutsinglelayer)
        continue;

      // check each panel to see if it has a bad hit pattern
      bool is_shower = false;
      int totalsamegap = 0;
      for (size_t i=0;i<straws.size();i++){
        if (is_shower)
          break;
        if (straws[i].size() > 0){
          std::sort(straws[i].begin(), straws[i].end());
          int layer = -1;
          bool crossed = false;
          for (size_t j=0;j<straws[i].size();j++){
            if (j == 0)
              layer = straws[i][j] % 2;
            else{
              int layer2 = straws[i][j] % 2;
              if (layer2 == layer){
                int gap = 0;
                for (size_t k=straws[i][j-1]+2;k<straws[i][j];k+=2){
                  StrawId sid(i/6,i%6,k);
                  if (!trackerStatus.noSignal(sid) && !trackerStatus.suppress(sid)) {
                    gap += 1;
                  }
                }
                totalsamegap += gap;
                if (gap > _maxSameGap || totalsamegap > _maxTotalSameGap){
                  is_shower = true;
                  break;
                }
              }else{
                // if already crossed, then
                if (crossed){
                  is_shower = true;
                  break;
                }else{
                  crossed = true;
                  int gap = 0;
                  for (size_t k=straws[i][j-1]+1;k<straws[i][j];k++){
                    StrawId sid(i/6,i%6,k);
                    if (!trackerStatus.noSignal(sid) && !trackerStatus.suppress(sid)) {
                      gap += 1;
                    }
                  }
                  if (gap > _maxCrossingGap){
                    is_shower = true;
                    break;
                  }
                }
              }
              layer = layer2;
            }
          }
        }
      }
      if (is_shower)
        continue;

      TimeClusterCollection*  tcol  = tccol_filter.get();
      tcol->push_back(tclust);
    }

    event.put(std::move(tccol_filter));
  }
}

using mu2e::CosmicShowerFilter;
DEFINE_ART_MODULE(CosmicShowerFilter)
