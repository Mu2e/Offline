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
#include "TTrackerGeom/inc/TTracker.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"

// art includes.
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

// From the art tool-chain
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C++ includes.
#include <iostream>

//#include "PerfLib/inc/perflib.hh"
//perf::PerfStats g_perf("FlagStrawHits new 200") ;

using namespace std;

namespace mu2e {

  class FlagStrawHits : public art::EDProducer {

  public:
    explicit FlagStrawHits(fhicl::ParameterSet const& pset);
    // Accept compiler written d'tor.

    void produce( art::Event& e);
    virtual void beginRun( art::Run& run );

  private:
   void initStrawToGlobalPanelMap();
      
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
    
    size_t _globalpanel_count;    
    vector<size_t> _map_strawindex_to_globalpanel; //straw to global panel map
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

void FlagStrawHits::initStrawToGlobalPanelMap(){
    const Tracker& tracker = getTrackerOrThrow();
    const TTracker& tt = dynamic_cast<const TTracker&>(tracker);

    size_t nplanes = tt.nPlanes();
    size_t npanels = tt.getPlane(0).nPanels();
    
    _globalpanel_count=nplanes*npanels;
    
    auto const&straws =tt.getAllStraws();
       
    _map_strawindex_to_globalpanel.reserve(straws.size());
    
    for(auto const&straw:straws){
      size_t iplane = straw.id().getPlane();
      size_t ipnl = straw.id().getPanel();
      size_t global_panel = ipnl + iplane*npanels;      
      _map_strawindex_to_globalpanel[straw.index().asInt()]=global_panel;
    }
}
void FlagStrawHits::beginRun( art::Run&){
  initStrawToGlobalPanelMap();
}


struct free_delete{
    void operator()(void* x) { free(x); }
};

  void
  FlagStrawHits::produce(art::Event& event) {
   //  g_perf.read_begin_counters_inlined();

    if ( _debug > 1 ) cout << "FlagStrawHits: produce() begin; event " << event.id().event() << endl;
      
    art::Handle<mu2e::StrawHitCollection> shcolH;    
    event.getByLabel(_shLabel,shcolH);
    
    const StrawHitCollection* shcol = shcolH.product();
    
    if(shcol == 0){
      throw cet::exception("RECO")<< "No StrawHit collection found for label " <<  _shLabel << endl;
    }

    if(shcol->size() >= INT_MAX)
	throw cet::exception("RECO")<<"mu2e::FlagStrawHits: shcol->size() >= INT_MAX" << endl;                 

    const StrawHitPositionCollection* shpcol(0);
    
    if(_shpLabel.length()!=0){      
      art::Handle<mu2e::StrawHitPositionCollection> shpcolH;      
      event.getByLabel(_shpLabel,shpcolH);      
      shpcol = shpcolH.product();      
      if(shpcol == 0){
	throw cet::exception("RECO") << "No StrawHitPosition collection found for label " <<  _shpLabel << endl;
      }
    }
  
    const size_t nsh = shcol->size();
    
    // create the output collection
    unique_ptr<StrawHitFlagCollection> shfcol(new StrawHitFlagCollection);
    shfcol->resize(nsh);

    // A more efficient algorithm is to first find all big hits, then
    // loop over only those in the double loop.  It avoids the search for indices FIXME!!
    vector<vector<int> > hits_by_panel(_globalpanel_count,vector<int>());

    for(auto & panel:hits_by_panel)
	panel.reserve(nsh>>6);
       
    for (size_t ish=0;ish<nsh;++ish){      
      StrawHit const& sh = shcol->at(ish);      
      const auto strawIndex=sh.strawIndex().asInt();      
      const auto accepted_ctE=sh.energyDep() > _ctE;

      StrawHitFlag& flag= shfcol->at(ish);
                 
      if(sh.energyDep() > _minE && sh.energyDep() < _maxE)
        flag.merge(StrawHitFlag::energysel);
      
      if(sh.time() > _minT && sh.time() < _maxT)
        flag.merge(StrawHitFlag::timesel);

      if(shpcol != 0){
       // merge with the position flag if that's present
	flag.merge(shpcol->at(ish).flag());	
        StrawHitPosition const& shp = shpcol->at(ish);	
        auto istereo = shp.flag().hasAllProperties(StrawHitFlag::stereo) ? 0 : 1;        
	auto rad = shp.pos().perp();
        
	if(rad > _minR && rad < _maxR[istereo])
          flag.merge(StrawHitFlag::radsel);
      }

      hits_by_panel[_map_strawindex_to_globalpanel[strawIndex]].push_back(accepted_ctE?ish:-ish);
  }
  {    
    std::unique_ptr<bool, free_delete>buffer((bool*)calloc(nsh<<1,sizeof(bool)));
    bool*const is_ct_straw=buffer.get();
 
    Tracker const& tracker = getTrackerOrThrow();
   
    for (auto const& panel:hits_by_panel){ // loop over panels
      for (auto ish : panel){
	if(ish<0) continue;
 
        StrawHit const& sh = shcol->at(ish);        
	Straw const& straw = tracker.getStraw( sh.strawIndex() );

          for (auto jsh : panel){
            if (ish == jsh) continue;	    
            	    
	    const size_t ush = jsh>0?jsh:-jsh;
	    
	    StrawHit const& sh2 = shcol->at(ush);
	    
	    if (sh2.time()-sh.time() > _ctMinT && sh2.time()-sh.time() < _ctMaxT){  
              if (!is_ct_straw[ush] && straw.isSamePreamp(sh2.strawIndex())){
		is_ct_straw[ush]=true;
		shfcol->at(ush).merge(StrawHitFlag::elecxtalk);
	      }
	      if (!is_ct_straw[ush+1] && straw.isNearestNeighbour(sh2.strawIndex())){
                is_ct_straw[ush+1]=true;
	        shfcol->at(ush).merge(StrawHitFlag::strawxtalk);
	      }
            }
            
          } // end loop over possible neighbors          
      } // end loop over hits in panel
    } // end loop over panels

  }
    event.put(move(shfcol));

 //  g_perf.read_end_counters_inlined();    

  } // end FlagStrawHits::produce.
  
} // end namespace mu2e

using mu2e::FlagStrawHits;
DEFINE_ART_MODULE(FlagStrawHits)

