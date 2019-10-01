//
//  Filter for selecting events with a mu-bunch intensity above a given threshold
//  Original author: G. Pezzullo
//
// framework
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "RecoDataProducts/inc/TrkFitFlag.hh"
#include "RecoDataProducts/inc/TriggerInfo.hh"
#include "fhiclcpp/ParameterSet.h"
#include "BFieldGeom/inc/BFieldManager.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
// data
#include "MCDataProducts/inc/ProtonBunchIntensity.hh"

using namespace CLHEP;
// c++
#include <string>
#include <vector>
#include <iostream>
#include <memory>

using namespace std;

namespace mu2e
{
  class BunchIntensityFilter : public art::EDFilter
  {
  public:
    explicit BunchIntensityFilter(fhicl::ParameterSet const& pset);
    virtual bool filter  (art::Event& event) override;
    virtual bool beginRun(art::Run&   run   );
    virtual bool endRun  ( art::Run& run ) override;

  private:
    art::InputTag _pbiTag;
    double        _minIntensity;
    int           _debug;
    // counters
    unsigned      _nevt, _npass;
  };

  BunchIntensityFilter::BunchIntensityFilter(fhicl::ParameterSet const& pset) :
    art::EDFilter{pset},
    _pbiTag       (pset.get<art::InputTag>("protonBunchIntensityTag","protonBunchIntensity")),
    _minIntensity (pset.get<double>       ("minIntensity", 3.5e7)),
    _debug        (pset.get<int>          ("debugLevel",0)),
    _nevt(0), _npass(0){}

  bool BunchIntensityFilter::beginRun(art::Run & run){
    return true;
  }

  bool BunchIntensityFilter::filter(art::Event& evt){
    // create output
    unique_ptr<TriggerInfo> triginfo(new TriggerInfo);
    ++_nevt;
    bool retval(false); // preset to fail
    // find the collection
    auto   pbiH = evt.getValidHandle<ProtonBunchIntensity>(_pbiTag);
    const  ProtonBunchIntensity* pbi = pbiH.product();

    if (pbi->intensity() >= _minIntensity) {
      retval = true;
      ++_npass;
    }

    return retval;
  }

  bool BunchIntensityFilter::endRun( art::Run& run ) {
    if(_debug > 0 && _nevt > 0){
      cout << moduleDescription().moduleLabel() << " paassed " <<  _npass << " events out of " << _nevt << " for a ratio of " << float(_npass)/float(_nevt) << endl;
    }
    return true;
  }
}
using mu2e::BunchIntensityFilter;
DEFINE_ART_MODULE(BunchIntensityFilter);
