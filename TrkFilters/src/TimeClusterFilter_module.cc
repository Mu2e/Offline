//
//  Filter for selecting good time cluster: this is part of the track trigger
//  Original author: Dave Brown (LBNL) 3/1/2017
//
// framework
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
// mu2e
// data
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/RecoDataProducts/inc/TriggerInfo.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
// c++
#include <iostream>
#include <memory>
#include <map>


namespace mu2e
{
  class TimeClusterFilter : public art::EDFilter
  {
    public:
      struct Config{
        using Name    = fhicl::Name;
        using Comment = fhicl::Comment;
        fhicl::Atom<art::InputTag>      timeClusterCollection{    Name("timeClusterCollection"),      Comment("TimeClusterCollection label") };
        fhicl::Atom<art::InputTag>      comboHitCollection   {    Name("ComboHitCollection"),         Comment("ComboHitCollection for timecluster"),"TTmakePH"};
        fhicl::Atom<bool>               requireCaloCluster   {    Name("requireCaloCluster"),         Comment("Require caloCluster") };
        fhicl::Atom<unsigned>           minNStrawHits        {    Name("minNStrawHits"),                   Comment("minNStrawHits")};
        fhicl::Atom<unsigned>           minNPlanes           {    Name("minNPlanes"),                 Comment("Minimum planes hit"), 0};
        fhicl::Atom<unsigned>           minPlaneSpan         {    Name("minPlaneSpan"),               Comment("Minimum difference between first and last hit plane"), 0};
        fhicl::Atom<int>                debugLevel           {    Name("debugLevel"),                 Comment("Debug"),0 };
        fhicl::Atom<bool>               noFilter             {    Name("noFilter"),                 Comment("Don't filter anything"),0 };
      };

      using Parameters = art::EDFilter::Table<Config>;

      explicit TimeClusterFilter(const Parameters& config);

    private:
      bool filter(art::Event& event) override;
      bool endRun(art::Run& run ) override;

      art::InputTag _tcTag;
      art::InputTag _chTag;
      bool          _hascc; // Calo Cluster
      unsigned      _minnhits;
      unsigned      _minnplanes;
      unsigned      _minplanespan;
      int           _debug;
      // counters
      unsigned      _nevt, _npass;
      bool          _noFilter;
  };

  TimeClusterFilter::TimeClusterFilter(const Parameters& conf)
    : art::EDFilter{conf},
    _tcTag   (conf().timeClusterCollection()),
    _chTag   (conf().comboHitCollection()),
    _hascc   (conf().requireCaloCluster()),
    _minnhits(conf().minNStrawHits()),
    _minnplanes(conf().minNPlanes()),
    _minplanespan(conf().minPlaneSpan()),
    _debug   (conf().debugLevel()),
    _nevt    (0),
    _npass   (0),
    _noFilter(conf().noFilter())
    {
      produces<TriggerInfo>();
    }

  bool TimeClusterFilter::filter(art::Event& evt){
    // create output
    std::unique_ptr<TriggerInfo> triginfo(new TriggerInfo);
    ++_nevt;
    bool retval(false); // preset to fail
    // find the collection
    auto tcH = evt.getValidHandle<TimeClusterCollection>(_tcTag);
    const TimeClusterCollection* tccol = tcH.product();
    auto chH = evt.getValidHandle<ComboHitCollection>(_chTag);
    const ComboHitCollection* chcol = chH.product();
    // loop over the collection: if any pass the selection, pass this event
    for(auto itc = tccol->begin();itc != tccol->end(); ++itc) {
      auto const& tc = *itc;
      if(_debug > 2){
        std::cout << moduleDescription().moduleLabel() << " nStrawHits = " << tc.nStrawHits() << " t0 = " << tc.t0().t0() << std::endl;
      }
      bool planeok = true;
      if (_minnplanes > 0 || _minplanespan > 0){
        planeok = false;
        int minplane = StrawId::_nplanes;
        int maxplane = 0;
        unsigned nplanes = 0;
        std::map<unsigned, unsigned> pmap;
        for (size_t ish = 0; ish < tc.hits().size(); ++ish){
          auto const& ch = chcol->at(tc.hits()[ish]);
          auto iplane = ch.strawId().plane();
          if (iplane < minplane) minplane = iplane;
          if (iplane > maxplane) maxplane = iplane;
          if (pmap.find(iplane) == pmap.end()){
            pmap[iplane] = 1;
            nplanes++;
          }
        }
        unsigned planespan = 0;
        if (maxplane > minplane) planespan = maxplane-minplane;
        if (nplanes >= _minnplanes && planespan >= _minplanespan){
          planeok = true;
        }
      }
      if( (!_hascc || tc.caloCluster().isNonnull()) &&
          tc.nStrawHits() >= _minnhits && planeok) {
        retval = true;
        ++_npass;
        // Fill the trigger info object
        // associate to the hit cluster which triggers.  Note there may be other hit clusters which also pass the filter
        // but filtering is by event!
        size_t index = std::distance(tccol->begin(),itc);
        triginfo->_hitClusters.push_back(art::Ptr<TimeCluster>(tcH,index));

        if(_debug > 1){
          std::cout << moduleDescription().moduleLabel() << " passed event " << evt.id() << std::endl;
        }
      }
    }
    evt.put(std::move(triginfo));
    if (!_noFilter){
      return retval;
    }else {
      return true;
    }
  }

  bool TimeClusterFilter::endRun( art::Run& run ) {
    if(_debug > 0 && _nevt > 0){
      std::cout << moduleDescription().moduleLabel() << " passed " << _npass << " events out of " << _nevt << " for a ratio of " << float(_npass)/float(_nevt) << std::endl;
    }
    return true;
  }
}
using mu2e::TimeClusterFilter;
DEFINE_ART_MODULE(TimeClusterFilter)
