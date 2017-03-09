//
//  Filter for selecting good time cluster: this is part of the track trigger
//  Original author: Dave Brown (LBNL) 3/1/2017
//
// framework
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
// mu2e
// data
#include "RecoDataProducts/inc/TimeCluster.hh"
// c++
#include <iostream>

using namespace std;

namespace mu2e 
{
  class TimeClusterFilter : public art::EDFilter 
  {
    public:
      explicit TimeClusterFilter(fhicl::ParameterSet const& pset);
      virtual bool filter(art::Event& event) override;
      virtual void endJob() override;

    private:
      art::InputTag _tcTag;
      bool _hascc; // Calo Cluster
      unsigned _minnhits;
      double _mintime, _maxtime;
      // counters
      unsigned _nevt, _npass;
  };

  TimeClusterFilter::TimeClusterFilter(fhicl::ParameterSet const& pset) :
    _tcTag(pset.get<art::InputTag>("TimeClusterCollection","TimeClusterFinder")),
    _hascc(pset.get<bool>("RequireCaloCluster",false)),
    _minnhits(pset.get<unsigned>("MinNHits",11)),
    _mintime(pset.get<double>("MinTime",500.0)),
    _maxtime(pset.get<double>("MaxTime",1695.0)) ,
    _nevt(0), _npass(0)
  {
  }

  bool TimeClusterFilter::filter(art::Event& evt){
    ++_nevt;
    bool retval(false); // preset to fail
// find the collection
    auto tcH = evt.getValidHandle<TimeClusterCollection>(_tcTag);
    const TimeClusterCollection* tccol = tcH.product();
// loop over the collection: if any pass the selection, pass this event
    for(auto const& tc : *tccol ) {
      if( (!_hascc || tc.caloCluster().isNonnull()) &&
	  tc.hits().size() >= _minnhits &&
	  tc.t0().t0() > _mintime && tc.t0().t0() < _maxtime) {
	retval = true;
	++_npass;
// should create a trigger product output FIXME
	break;
      }
    }
    return retval;
  }

  void TimeClusterFilter::endJob() {
    if(_nevt > 0){
      cout << "TimeClusterFilter passed " <<  _npass << " events out of " << _nevt << " for a ratio of " << float(_npass)/float(_nevt) << endl;
    }

  }
}
