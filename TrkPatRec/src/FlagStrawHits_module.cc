//
// A module to flag StrawHits for track reconstruction and delta ray 
// identification
//
// $Id: FlagStrawHits_module.cc,v 1.13 2014/05/28 20:29:28 brownd Exp $
// $Author: brownd $
// $Date: 2014/05/28 20:29:28 $
// 
//  Original Author: David Brown, LBNL
//  

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
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

  class FlagStrawHits : public art::EDProducer {

  public:
    explicit FlagStrawHits(fhicl::ParameterSet const& pset);
    // Accept compiler written d'tor.

    void produce( art::Event& e);

  private:

    // Diagnostics level.
    int _diag;
    int _debug;
    // input collection labels
    string _shLabel;
    string _shpLabel; 
    // Parameters; these are vectors for the number of quality selections
    double _minE, _maxE;  // minimum and maximum charge (units??)
    double _minT, _maxT; // minimum and maximum hit time
    double _minR;
    vector<double> _maxR; // minimum and maximum transverse radius
  };

  FlagStrawHits::FlagStrawHits(fhicl::ParameterSet const& pset) :
    _diag(pset.get<int>("diagLevel",0)),
    _debug(pset.get<int>("debugLevel",0)),
    _shLabel(pset.get<string>("StrawHitCollectionLabel","makeSH")),
    _shpLabel(pset.get<string>("StrawHitPositionCollectionLabel","MakeStereoHits")),
    _minE(pset.get<double>("minimumEnergy",0.0)),
    _maxE(pset.get<double>("maximumEnergy",0.003)),
    _minT(pset.get<double>("minimumTime",500)),
    _maxT(pset.get<double>("maximumTime",2000)),
    _minR(pset.get<double>("minimumRadius",395.0)),
    _maxR(pset.get<vector<double> >("maximumRadius",vector<double>{650,650}))
    {
      if(_maxR.size() != 2)
	throw cet::exception("RECO")<<"mu2e::FlagStrawHits: illegal maximumRadius size specified" << endl;
      produces<StrawHitFlagCollection>();
    }
  void
  FlagStrawHits::produce(art::Event& event) {

    if ( _debug > 1 ) cout << "FlagStrawHits: produce() begin; event " << event.id().event() << endl;

    art::Handle<mu2e::StrawHitCollection> shcolH;
    if(event.getByLabel(_shLabel,shcolH));
    const StrawHitCollection* shcol = shcolH.product();
    if(shcol == 0){
      if(_debug > 0) cout << "No StrawHit collection found for label " <<  _shLabel << endl;
      return;
    }
    const StrawHitPositionCollection* shpcol(0);
    if(_shpLabel.length() >0){
      art::Handle<mu2e::StrawHitPositionCollection> shpcolH;
      if(event.getByLabel(_shpLabel,shpcolH));
      shpcol = shpcolH.product();
      if(shpcol == 0){
	if(_debug > 0) cout << "No StrawHitPosition collection found for label " <<  _shpLabel << endl;
	return;
      }
    }
  
    size_t nsh = shcol->size();
    // create the output collection
    unique_ptr<StrawHitFlagCollection> shfcol(new StrawHitFlagCollection);
    shfcol->reserve(nsh);
    for(size_t ish=0;ish<nsh;++ish){
      StrawHit const& sh = shcol->at(ish);
      StrawHitFlag flag;
      if(sh.energyDep() > _minE && sh.energyDep() < _maxE)
	flag.merge(StrawHitFlag::energysel);
      if(sh.time() > _minT && sh.time() < _maxT)
	flag.merge(StrawHitFlag::timesel);
      if(shpcol != 0){
	StrawHitPosition const& shp = shpcol->at(ish);
	unsigned istereo = shp.flag().hasAllProperties(StrawHitFlag::stereo) ? 0 : 1;
	double rad = shp.pos().perp();
	if(rad > _minR && rad < _maxR[istereo])
	  flag.merge(StrawHitFlag::radsel);
      }
      shfcol->push_back(flag);
    }
    event.put(move(shfcol));

  } // end FlagStrawHits::produce.
} // end namespace mu2e

using mu2e::FlagStrawHits;
DEFINE_ART_MODULE(FlagStrawHits)

