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
#include "RecoDataProducts/inc/StrawHitCollection.hh"
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
    double _ctE; // minimum charge to flag neighbors as cross talk
    double _ctMinT, _ctMaxT; // time relative to proton hit to flag cross talk (ns)
    double _minT, _maxT; // minimum and maximum hit time
    double _minR;
    vector<double> _maxR; // minimum and maximum transverse radius
  };

  FlagStrawHits::FlagStrawHits(fhicl::ParameterSet const& pset) :
    _diag(pset.get<int>("diagLevel",0)),
    _debug(pset.get<int>("debugLevel",0)),
    _shLabel(pset.get<string>("StrawHitCollectionLabel","makeSH")),
    _shpLabel(pset.get<string>("StrawHitPositionCollectionLabel","MakeStereoHits")),
    _minE(pset.get<double>("minimumEnergy",0.0)), // MeV
    _maxE(pset.get<double>("maximumEnergy",0.0035)), // MeV
    _ctE(pset.get<double>("crossTalkEnergy",0.007)), // MeV
    _ctMinT(pset.get<double>("crossTalkMinimumTime",-1)), // nsec
    _ctMaxT(pset.get<double>("crossTalkMaximumTime",100)), // nsec
    _minT(pset.get<double>("minimumTime",500)), // nsec
    _maxT(pset.get<double>("maximumTime",2000)), // nsec
    _minR(pset.get<double>("minimumRadius",395.0)), // mm
    _maxR(pset.get<vector<double> >("maximumRadius",vector<double>{650,650})) // mm
    {
      if(_maxR.size() != 2)
	throw cet::exception("RECO")<<"mu2e::FlagStrawHits: illegal maximumRadius size specified" << endl;
      produces<StrawHitFlagCollection>();
    }

  void
  FlagStrawHits::produce(art::Event& event) {
    const Tracker& tracker = getTrackerOrThrow();

    if ( _debug > 1 ) cout << "FlagStrawHits: produce() begin; event " << event.id().event() << endl;

    art::Handle<mu2e::StrawHitCollection> shcolH;
    if(event.getByLabel(_shLabel,shcolH));
    const StrawHitCollection* shcol = shcolH.product();
    if(shcol == 0){
      throw cet::exception("RECO")<< "No StrawHit collection found for label " <<  _shLabel << endl;
    }
    const StrawHitPositionCollection* shpcol(0);
    if(_shpLabel.length() >0){
      art::Handle<mu2e::StrawHitPositionCollection> shpcolH;
      if(event.getByLabel(_shpLabel,shpcolH));
      shpcol = shpcolH.product();
      if(shpcol == 0){
	throw cet::exception("RECO") << "No StrawHitPosition collection found for label " <<  _shpLabel << endl;
      }
    }
  
    size_t nsh = shcol->size();
    // create the output collection
    unique_ptr<StrawHitFlagCollection> shfcol(new StrawHitFlagCollection);
    shfcol->reserve(nsh);
// A more efficient algorithm is to first find all big hits, then
// loop over only those in the double loop.  It avoids the search for indices FIXME!!

    vector<StrawIndex> ct_straws_neighbor;
    vector<StrawIndex> ct_straws_preamp;
    
    for (size_t ish=0;ish<nsh;++ish){
      StrawHit const& sh = shcol->at(ish);
      const Straw& straw = tracker.getStraw( sh.strawIndex() );
      if (sh.energyDep() >= _ctE){
        for (size_t jsh=0;jsh<nsh;++jsh){
          StrawHit const& sh2 = shcol->at(jsh);
          if (sh2.time()-sh.time() > _ctMinT && sh2.time()-sh.time() < _ctMaxT){
            if (straw.isSamePreamp(sh2.strawIndex()))
              ct_straws_preamp.push_back(sh2.strawIndex());
            if (straw.isNearestNeighbour(sh2.strawIndex()))
              ct_straws_neighbor.push_back(sh2.strawIndex());
          }
        }
      }
    }

    for(size_t ish=0;ish<nsh;++ish){
      StrawHit const& sh = shcol->at(ish);
      StrawHitFlag flag;
      if(sh.energyDep() > _minE && sh.energyDep() < _maxE)
        flag.merge(StrawHitFlag::energysel);

      if (find(ct_straws_neighbor.begin(),ct_straws_neighbor.end(),sh.strawIndex()) != ct_straws_neighbor.end())
        flag.merge(StrawHitFlag::strawxtalk);
      if (find(ct_straws_preamp.begin(),ct_straws_preamp.end(),sh.strawIndex()) != ct_straws_preamp.end())
        flag.merge(StrawHitFlag::elecxtalk);

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

