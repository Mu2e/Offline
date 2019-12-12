//
//  Filter for selecting good seed (chisq) track fits: this is part of the track trigger
//  Original author: Dave Brown (LBNL) 3/1/2017
//
// framework
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
//Author: S Middleton
//Date: Dec 2019
//Purpose: For Devloping Straight Track Trigger
#include "RecoDataProducts/inc/TrkFitFlag.hh"
#include "fhiclcpp/ParameterSet.h"

#include "RecoDataProducts/inc/CosmicTrackSeed.hh"
#include "RecoDataProducts/inc/TriggerInfo.hh"

using namespace CLHEP;
// c++
#include <string>
#include <vector>
#include <iostream>
#include <memory>

using namespace std;

namespace mu2e
{
  class SeedFilter : public art::EDFilter
  {
  public:
    explicit SeedFilter(fhicl::ParameterSet const& pset);
    virtual bool filter(art::Event& event) override;
    virtual bool endRun( art::Run& run ) override;

  private:

    art::InputTag   _cosmicTag;

    int _minnsh;
    int _minnch;    
    int _minNHitsTimeCluster;
/*
    unsigned n_outliers;
    unsigned maxniter;

*/
    TrkFitFlag      _goodcosmic; 
    TrkFitFlag      _convergedcosmic;
    std::string     _trigPath;
    int             _debug;
    // counters
    unsigned        _nevt, _npass;
  };

  SeedFilter::SeedFilter(fhicl::ParameterSet const& pset) :
    art::EDFilter{pset},
    _cosmicTag     (pset.get<art::InputTag>("CosmicTrackSeedCollection","KSFDeM")),
    _goodcosmic(pset.get<vector<string> >("comsicseedFitFlag",vector<string>{"HelixOK"})),
    _convergedcosmic(pset.get<vector<string> > ("comsicseedFitFlag", vector<string>{"HelixConverged"})),
    _minnsh (pset.get<int>   ("minnsh",8)),
    _minnch (pset.get<int>   ("minnch",8)),
    _trigPath  (pset.get<std::string>("triggerPath")),
    _debug     (pset.get<int>   ("debugLevel",0)),
    _nevt(0), _npass(0)
  {
    produces<TriggerInfo>();
  }

  bool SeedFilter::filter(art::Event& evt){
    unique_ptr<TriggerInfo> triginfo(new TriggerInfo);
    ++_nevt;
    bool   retval(false);
    // find the collection
    auto cosH = evt.getValidHandle<CosmicTrackSeedCollection>(_cosmicTag);
    const CosmicTrackSeedCollection* coscol = cosH.product();
    // loop over the collection: if any pass the selection, pass this event
    for(auto icos = coscol->begin(); icos != coscol->end(); ++icos) {
      auto const& cosmic = *icos;
     
      if( cosmic.status().hasAllProperties(_goodcosmic) && cosmic.status().hasAllProperties(_convergedcosmic) && cosmic.hits().size()>_minnch && cosmic.trkstrawhits().size() > _minnsh &&cosmic.timeCluster().isNonnull()){ 
       
        ++_npass;
        // Fill the trigger info object
        triginfo->_triggerBits.merge(TriggerFlag::track); //use the old track flag
        triginfo->_triggerPath = _trigPath;
        // associate to the helix which triggers.  Note there may be other helices which also pass the filter
        // but filtering is by event!
        size_t index = std::distance(coscol->begin(),icos);
        triginfo->_track = art::Ptr<CosmicTrackSeed>(cosH,index);
        if(_debug > 1){
          cout << moduleDescription().moduleLabel() << " passed event " << evt.id() << endl;
        }
        break;
      }
    }
    evt.put(std::move(triginfo));
    return retval;
  }

  bool CosmicSeedFilter::endRun( art::Run& run ) {
    if(_debug > 0 && _nevt > 0){
      cout << moduleDescription().moduleLabel() << " passed " <<  _npass << " events out of " << _nevt << " for a ratio of " << float(_npass)/float(_nevt) << endl;
    }
    return true;
  }
}
using mu2e::CosmicSeedFilter;
DEFINE_ART_MODULE(CosmicSeedFilter);
