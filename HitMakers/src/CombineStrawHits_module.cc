//
// A module to combine straw hits in the same panel.  This improves resolution and
// reduces combinatorix for downstream modules.
//  Original Author: David Brown, LBNL
//  

// Mu2e includes.
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/StrawHitPosition.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
// art includes.
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

// From the art tool-chain
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
// boost
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/weighted_variance.hpp> 

// C++ includes.
#include <iostream>
#include <float.h>

using namespace std;

namespace mu2e {

  // struct for combining hits
  struct CHInfo {
    StrawHitIndex _index; // index to straw hit
    uint16_t _straw; // straw number
    float _wd; // distance along wire from center
    float _werr2; // error for this distance (squared)
    float _time; // hit time
  };

  class CombineStrawHits : public art::EDProducer {

    public:
      explicit CombineStrawHits(fhicl::ParameterSet const& pset);

      void produce( art::Event& e);

    private:

      // configuration
      int _debug;
      // event object Tags
      art::InputTag   _shTag, _shpTag, _shfTag;
      // input event collection cache
      const StrawHitPositionCollection*     _shpcol;
      // Parameters
      StrawHitFlag _shsel; // flag selection
      StrawHitFlag _shmask; // flag anti-selection 
      double _maxdt; // maximum time separation between hits
      double _maxwdchi; // maximum wire distance separation chi
      int _maxds; // maximum straw number difference
  };

  CombineStrawHits::CombineStrawHits(fhicl::ParameterSet const& pset) :
    // Parameters
    _debug(pset.get<int>("debugLevel",0)),
    _shTag		(pset.get<art::InputTag>("StrawHitCollection","makeSH")),
    _shpTag		(pset.get<art::InputTag>("StrawHitPositionCollection","MakeStrawHitPositions")),
    _shfTag	       (pset.get<art::InputTag>("StrawHitFlagCollection","FSHPreStereo")),
    _shsel(pset.get<vector<string> >("StrawHitSelectionBits",vector<string>{"EnergySelection","TimeSelection"} )),
    _shmask(pset.get<vector<string> >("StrawHitMaskBits",vector<string>{} )),
    _maxdt(pset.get<double>("MaxDt",40.0)), // nsec
    _maxwdchi(pset.get<double>("MaxWireDistDiffPull",4.0)), //units of resolution sigma
    _maxds(pset.get<int>("MaxDS",3)) // how far away 2 straws can be, in 0-95 numbering (including layers!!)
  {
    produces<ComboHitCollection>();
  }

  void CombineStrawHits::produce(art::Event& event) {
    // Get a reference to T trackers
    const Tracker& tracker = getTrackerOrThrow();
    const TTracker& tt = dynamic_cast<const TTracker&>(tracker);

    // find event data
    if( !findData(event) ){
      throw cet::exception("RECO")<<"mu2e::CombineStrawHits: collection missing or sizes don't match " << endl; 
    }
   // create output 
    unique_ptr<ComboHitCollection> chcol(new ComboHitCollection());
 
 // sort hits by panel
    std::array<std::vector<CHInfo>,240> panelhits;

    size_t nsh = _shcol->size();
    for(uint16_t ish=0;ish<nsh;++ish){
      // merge the flags
      StrawHitFlag shf = _shfcol->at(ish);
      shf.merge(_shpcol->at(ish).flag());
    // select hits based on flag
      if(shf.hasAllProperties(_shsel) && (!shf.hasAnyProperty(_shmask)) ){
	StrawHit const& hit = (*_shcol)[ish];
	StrawHitPosition const& hitpos = (*_shpcol)[ish];
	Straw const& straw = tt.getStraw(hit.strawIndex());
        int plane = straw.id().getPlane();
        int panel = straw.id().getPanel();
	unsigned index = plane*40+panel; // crude indexing, should be provided by StrawId, FIXME!
	SHinfo info;
	info._index = ish;
	info._straw = straw.id().getStraw();
	info._wd = hitpos.wireDist();
	info._werr2 = hitpos._wres*hitpos._wres;
	info._time = std::min(hit.time(StrawEnd::cal),hit.time(StrawEnd::hv));
	panelhits[index].push_back(info);
      }
    }
    // keep track of which hits are used as part of a combo hit
    std::array<bool,nsh> used{false};
    // loop over panels
    for(auto const& phits : panelhits ) {
    // loop over hit pairs in this panel
      for(size_t ihit=0;ihit < phits.size(); ++ihit){
	CHInfo const& ish = phits[ihit];
	if(!used[ish._index]){
	  // create a combo hit for every hit; initialize it with this hit
	  ComboHit combohit;
	  combohit._sh[combohit._nsh]=ish._index;
	  ++combohit._nsh;
	  used[ish._index] = true;
	  for(size_t jhit=ihit+1;jhit < phits.size(); ++jhit){
	    CHInfo const& jsh = phits[jhit];
	    if(!used[jsh._index]){
	      // require straws be near each other
	      int ds = abs( (int)ish._straw-(int)jsh._straw);
	      if(ds > 0 && ds <= _maxds ){
		// require times be consistent
		float dt = fabs(ish._time - jsh._time);
		if(dt < _maxdt){
		  // compute the chi of the differnce in wire positions
		  float wderr = sqrt(ish._werr2 + jsh._werr2);
		  float wdchi = fabs(ish._wd - jsh._wd)/wderr;
		  if(wdchi < _maxwdchi){
		  // add a neural net selection here someday for Offline use  FIXME!
		    // these hits match: add the 2nd to the combo hit
		    combohit._sh[combohit._nsh] = jsh._index;
		    ++combohit._nsh;
		    used[jsh._index] = true;
		  } // positions along wire
		}// consistent times
	      }// straw proximity 
	    } // 2nd hit not used
	  } // 2nd panel hit
	} // 1st hit not used
      } // 1st panel hit
      // compute floating point info for this combo hit and save it
      combineHits(combohit);
      chcol->push_back(combohit);
    } // panels
    // store data in the event
    event.put(std::move(chcol));
  }

  bool CombineStrawHits::findData(art::Event& evt){
    _shcol = 0; _shpcol = 0;  _shfcol = 0;
    auto shH = evt.getValidHandle<StrawHitCollection>(_shTag);
    _shcol = shH.product();
    auto shpH = evt.getValidHandle<StrawHitPositionCollection>(_shpTag);
    _shpcol = shpH.product();
    auto shfH = evt.getValidHandle<StrawHitFlagCollection>(_shfTag);
    _shfcol = shfH.product();
    return _shcol != 0 && _shpcol != 0 && _shpcol->size() == _shcol->size();
  }

// compute the combined properties of this combined hit
  void CombineStrawHits::combineHits(ComboHit& combohit) {
    combohit._wdir = (*_shpcol)[combohit._sh[0]].wdir();
    // if there's only 1 hit, take the info from the orginal collections
    if(combohit._nsh == 1){
      combohit._pos = (*_shpcol)[combohit._sh[0]].pos(); 
      combohit._wres = (*_shpcol)[combohit._sh[0]]._wres;
      combohit._tres = (*_shpcol)[combohit._sh[0]]._tres;
      combohit._time =  std::min((*_shcol)[combohit._sh[0]]._time[0], (*_shcol)[combohit._sh[0]]._time[1]);
      combohit._edep = (*_shcol)[combohit._sh[0]].energyDep();
    } else {

      accumulator_set<float, stats<tag::mean>> tacc;
      for(unsigned ish = 0; ish < combohit._nsh; ++ish){
	// get back the original information
	StrawHit const& sh = (*_shcol)[combohit._sh[ish]];
	StrawHitPosition const& shp = (*_shpcol)[combohit._sh[ish]];
	tacc(std::min(sh._tcal,sh._thv));
	wacc(shp._wdist);
      }
    }
  }
} // end namespace mu2e

  using mu2e::CombineStrawHits;
DEFINE_ART_MODULE(CombineStrawHits)

