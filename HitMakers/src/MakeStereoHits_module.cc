//
// A module to create simple stereo hits out of StrawHits.  This can work
// with either tracker.  StrawHit selection is done by flagging in an upstream module
//
// $Id: MakeStereoHits_module.cc,v 1.1 2013/01/26 18:19:17 brownd Exp $
// $Author: brownd $
// $Date: 2013/01/26 18:19:17 $
// 
//  Original Author: David Brown, LBNL
//  

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StereoHitCollection.hh"
#include "RecoDataProducts/inc/DFlagCollection.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"

// art includes.
#include "art/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

// From the art tool-chain
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C++ includes.
#include <iostream>

using namespace std;

namespace mu2e {

  class MakeStereoHits : public art::EDProducer {

  public:
    explicit MakeStereoHits(fhicl::ParameterSet const& pset);
    // Accept compiler written d'tor.

    void produce( art::Event& e);

  private:

    // Diagnostics level.
    int _diagLevel;
    // Name of the StrawHit collection
    std::string _strawHitsLabel;
    // name of corresponding DFlag collection
    std::string _strawHitFlagsLabel;
  // Parameters
    int _hqual;
    unsigned _hprop;
    double _maxDz; // maximum z sepration
    double _maxDt; // maximum time separation between hits
    double _maxDE; // maximum deposited energy deference: this excludes inconsistent hits
    double _maxDot; // maximum dot product magnitude between wire directions
    bool _allPairs; // whether to make all possible stereo pairs, or only use each hit once
    double _maxTDNsig; // maximum # of TimeDivision sigmas past the active edge of a wire to allow making stereo hits
// local cache of event data
    const StrawHitCollection* _strawhits; 
    const DFlagCollection* _shflags;
// selection function
    bool select(size_t ish) const;
  };

  MakeStereoHits::MakeStereoHits(fhicl::ParameterSet const& pset) :

    // Parameters
    _diagLevel(pset.get<int>("diagLevel",0)),
    _strawHitsLabel(pset.get<string>("strawHitsLabel","makeSH")),
    _strawHitFlagsLabel(pset.get<string>("strawHitFlagsLabel","FlagStrawHits")),
    _hqual(pset.get<int>("HitQuality",1)),
    _hprop(pset.get<unsigned>("HitProperty",0)),
    _maxDz(pset.get<double>("maxDz",35.0)), // mm
    _maxDt(pset.get<double>("maxDt",60.0)), // nsec
    _maxDE(pset.get<double>("maxDE",0.1)), // dimensionless
    _maxDot(pset.get<double>("maxDot",0.9)), // must be < 1.0
    _allPairs(pset.get<bool>("allPairs",true)), // make all possible pairs
    _maxTDNsig(pset.get<double>("maxTDNsig",3.0)),
    _strawhits(0), _shflags(0){

    // Tell the framework what we make.
    produces<StereoHitCollection>();
  }

  
  void
  MakeStereoHits::produce(art::Event& event) {

    if ( _diagLevel > 1 ) cout << "MakeStereoHits: produce() begin; event " << event.id().event() << endl;

    // Get a reference to one of the L or T trackers.
    // Throw exception if not successful.
    const Tracker& tracker = getTrackerOrThrow();

   // Handle to the conditions service
    ConditionsHandle<TrackerCalibrations> tcal("ignored");

    art::Handle<mu2e::StrawHitCollection> strawhitsH;
    if(event.getByLabel(_strawHitsLabel,strawhitsH))
    _strawhits = strawhitsH.product();
    if(_strawhits == 0){
      if(_diagLevel > 0) cout << "No StrawHit collection found for label " <<  _strawHitsLabel << endl;
      return;
    }
    art::Handle<mu2e::DFlagCollection> shflagsH;
    if(event.getByLabel(_strawHitFlagsLabel,shflagsH))
    _shflags = shflagsH.product();

    if(_shflags != 0){
      if(_strawhits->size() != _shflags->size())
	throw cet::exception("RECO")<<"mu2e::MakeStereoHits: inconsistent input collections" << endl;
    } else if (_diagLevel > 0 && _strawHitFlagsLabel != "None"){
      cout << "No StrawHit flag collection found for label " << _strawHitFlagsLabel << endl;
    }
    // select the straws, and save TD info
    std::vector<size_t> selsh;
    selsh.reserve(_strawhits->size());
    std::vector<SHInfo> shinfos;
    shinfos.reserve(_strawhits->size());
    for(size_t ish=0;ish<_strawhits->size();++ish){
      if(select(ish)){
	selsh.push_back(ish);
	SHInfo shinfo;
	Straw const& straw= tracker.getStraw((*_strawhits)[ish].strawIndex());
	tcal->StrawHitInfo(straw,(*_strawhits)[ish],shinfo);
	shinfos.push_back(shinfo);
      }
    }
    size_t nsh = selsh.size();
    // create the output collection
    auto_ptr<StereoHitCollection> stereoHits(new StereoHitCollection);
    stereoHits->reserve(1.5*nsh);
    // double loop over selected straw hits
    for(size_t ish=0;ish<nsh;++ish){
      StrawHit const& sh1 = (*_strawhits)[selsh[ish]];
      const Straw& straw1 = tracker.getStraw(sh1.strawIndex());
      for(size_t jsh=ish+1;jsh<nsh;++jsh){
	StrawHit const& sh2 = (*_strawhits)[selsh[jsh]];
	const Straw& straw2 = tracker.getStraw(sh2.strawIndex());
	if( fabs(sh1.time()-sh2.time()) < _maxDt // hits are close in time
	    && fabs(straw1.getMidPoint().z()-straw2.getMidPoint().z()) < _maxDz // z separation
	    && fabs(sh1.energyDep() - sh2.energyDep()) < _maxDE // deposited energy is roughly consistent (should compare dE/dx but have no path yet!)
	    && fabs(straw1.getDirection().dot(straw2.getDirection())) < _maxDot) // wires aren't paralell
	{
  // tentative stereo hit: this solves for the POCA
	  StereoHit sth(*_strawhits,tracker,selsh[ish],selsh[jsh]);
	  if( straw1.getDetail().activeHalfLength()-fabs(sth.wdist1()) > -straw1.getDetail().outerRadius() // stereo point is inside the active length
	      && straw2.getDetail().activeHalfLength()-fabs(sth.wdist2()) > -straw2.getDetail().outerRadius() // stereo point is inside the active length
	      && fabs((shinfos[ish]._pos-sth.pos()).dot(straw1.getDirection()))/shinfos[ish]._tdres<_maxTDNsig // stereo point is consistent with TD
	      && fabs((shinfos[jsh]._pos-sth.pos()).dot(straw2.getDirection()))/shinfos[jsh]._tdres<_maxTDNsig)
	  {
	    stereoHits->push_back(sth);
	  }
	}
      }
    }
    event.put(stereoHits);

  } // end MakeStereoHits::produce.

  bool 
  MakeStereoHits::select(size_t ish) const {
    return _shflags == 0 || (*_shflags)[ish].select(_hprop,_hqual);
  }

} // end namespace mu2e

using mu2e::MakeStereoHits;
DEFINE_ART_MODULE(MakeStereoHits)

