//
// A module to create simple stereo hits out of StrawHits.  This can work
// with either tracker.  StrawHit selection is done by flagging in an upstream module
//
// $Id: MakeStereoHits_module.cc,v 1.2 2013/03/08 04:30:23 brownd Exp $
// $Author: brownd $
// $Date: 2013/03/08 04:30:23 $
// 
//  Original Author: David Brown, LBNL
//  

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StereoHitCollection.hh"
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
    std::string _shLabel;
  // Parameters
    double _maxDz; // maximum z sepration
    double _maxDt; // maximum time separation between hits
    double _maxDE; // maximum deposited energy deference: this excludes inconsistent hits
    double _maxDot; // maximum dot product magnitude between wire directions
    bool _allPairs; // whether to make all possible stereo pairs, or only use each hit once
    double _maxTDNsig; // maximum # of TimeDivision sigmas past the active edge of a wire to allow making stereo hits
    double _maxChisq; // maximum # of TimeDivision consistency chisquared to allow making stereo hits
  };

  MakeStereoHits::MakeStereoHits(fhicl::ParameterSet const& pset) :

    // Parameters
    _diagLevel(pset.get<int>("diagLevel",0)),
    _shLabel(pset.get<string>("StrawHitCollectionLabel","makeSH")),
    _maxDz(pset.get<double>("maxDz",75.0)), // mm
    _maxDt(pset.get<double>("maxDt",60.0)), // nsec
    _maxDE(pset.get<double>("maxDE",0.95)), // dimensionless, from 0 to 1
    _maxDot(pset.get<double>("maxDot",0.98)), // must be < 1.0
    _allPairs(pset.get<bool>("allPairs",true)), // make all possible pairs
    _maxTDNsig(pset.get<double>("maxTDNsig",4.0)),
    _maxChisq(pset.get<double>("maxChisquared",9.0)){
    // Tell the framework what we make.
    produces<StereoHitCollection>();
    produces<StrawHitPositionCollection>();
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
    const StrawHitCollection* strawhits(0);
    if(event.getByLabel(_shLabel,strawhitsH))
      strawhits = strawhitsH.product();
    if(strawhits == 0){
      throw cet::exception("RECO")<<"mu2e::MakeStereoHits: No StrawHit collection found for label " <<  _shLabel << endl;
    }
    // create a collection of StrawHitPosition, and intialize them using the time division
    size_t nsh = strawhits->size();
    auto_ptr<StrawHitPositionCollection> shpos(new StrawHitPositionCollection);
    shpos->reserve(2*nsh);
    for(size_t ish=0;ish<nsh;++ish){
      shpos->push_back(StrawHitPosition(strawhits->at(ish),tracker,tcal));
    }
    // create the stereo hits
    auto_ptr<StereoHitCollection> stereohits(new StereoHitCollection);
    stereohits->reserve(1.5*nsh);
    // double loop over selected straw hits
    for(size_t ish=0;ish<nsh;++ish){
      StrawHit const& sh1 = strawhits->at(ish);
      const Straw& straw1 = tracker.getStraw(sh1.strawIndex());
      for(size_t jsh=ish+1;jsh<nsh;++jsh){
	StrawHit const& sh2 = strawhits->at(jsh);
	const Straw& straw2 = tracker.getStraw(sh2.strawIndex());
	if( fabs(sh1.time()-sh2.time()) < _maxDt // hits are close in time
	    && fabs(straw1.getMidPoint().z()-straw2.getMidPoint().z()) < _maxDz // z separation
	    && fabs(sh1.energyDep() - sh2.energyDep())/(sh1.energyDep()+sh2.energyDep()) < _maxDE // deposited energy is roughly consistent (should compare dE/dx but have no path yet!)
	    && fabs(straw1.getDirection().dot(straw2.getDirection())) < _maxDot) // wires aren't paralell
	{
  // tentative stereo hit: this solves for the POCA
	  StereoHit sth(*strawhits,tracker,ish,jsh);
	  if( straw1.getDetail().activeHalfLength()-fabs(sth.wdist1()) > -straw1.getDetail().outerRadius() // stereo point is inside the active length
	      && straw2.getDetail().activeHalfLength()-fabs(sth.wdist2()) > -straw2.getDetail().outerRadius()) { // stereo point is inside the active length
	    // compute difference between stereo points and TD prediction
	    double chi1 = ((*shpos)[ish].pos()-sth.pos()).dot(straw1.getDirection())/(*shpos)[ish].posRes(StrawHitPosition::phi);	      
	    double chi2 = ((*shpos)[jsh].pos()-sth.pos()).dot(straw2.getDirection())/(*shpos)[jsh].posRes(StrawHitPosition::phi);
	    if(fabs(chi1) <_maxTDNsig && fabs(chi2) < _maxTDNsig)
	    {
	      // compute chisquared
	      double chisq = chi1*chi1+chi2*chi2; 
	      if(chisq < _maxChisq){
		sth.setChisquared(chisq);
		stereohits->push_back(sth);
	      }
	    }
	  }
	}
      }
    }
// now, resolve the stereo hits to find the best position for each hit that particpates.  The algorithm is:
// 1) if a hit pairs with another hit in the same station, take the pair with the minimum plane separation
// 2) if a hit pairs with another hit in the same plane, take the pair with the minimum chisquared
    size_t nst = stereohits->size();
    for(size_t ish=0; ish<nsh;++ish){
      StrawHitPosition& shp = shpos->at(ish);
      double minchisq(100.0);
      SectorId::isep minsep(SectorId::apart);
      for(size_t ist=0;ist<nst;++ist){
	StereoHit const& sth = stereohits->at(ist);
	if(sth.hitIndex1() == ish || sth.hitIndex2() == ish) {
	  if(sth.sectorSeparation() < minsep ||
	      (sth.sectorSeparation() == minsep && sth.chisq() < minchisq)){
	    minsep = sth.sectorSeparation();
	    minchisq = sth.chisq();
	    shp = StrawHitPosition(sth);
	  }
	}
      }
    }
    event.put(stereohits);
    event.put(shpos);
  } // end MakeStereoHits::produce.

} // end namespace mu2e

using mu2e::MakeStereoHits;
DEFINE_ART_MODULE(MakeStereoHits)

