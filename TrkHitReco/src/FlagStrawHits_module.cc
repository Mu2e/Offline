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
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "RecoDataProducts/inc/StrawHitPosition.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
#include "TrkReco/inc/TrkTimeCalculator.hh"
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

using namespace std;

namespace mu2e {
  //  typedef pair<float,float> TimeRange;
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
    art::InputTag _shTag, _shpTag, _ccTag;
    const StrawHitCollection*             _shcol;
    const StrawHitPositionCollection*     _shpcol;
    const CaloClusterCollection*	  _cccol;
    bool  _usecc;
    bool   _useccxy;
    float _ccmine; // calo cluster energy minimum
    int    _ccminsize;// calo cluster minimum size given in crystal units
    float _ccmindt, _ccmaxdt;
    float _pitchangle;
    float _beta;
    float _dtoffset;//time offset between the strawHits and CaloCrystalHits
// Parameters; these are vectors for the number of quality selections
    float _minE, _maxE;  // minimum and maximum charge (units??)
    float _ctE; // minimum charge to flag neighbors as cross talk
    float _ctMinT, _ctMaxT; // time relative to proton hit to flag cross talk (ns)
    float _minT, _maxT; // minimum and maximum hit time
    float _minR, _maxR; // minimum and maximum transverse radius
    // T0 calcuator
    TrkTimeCalculator _ttcalc;
    float _hmaxtdrift; // cache of 1/2 the maximum drift time
  // helper functions
    bool goodCaloCluster(CaloCluster const& cc) const;
    //    bool goodTime(float hittime,vector<TimeRange> const& timeranges) const;
    bool findData(const art::Event& evt);
  };

  FlagStrawHits::FlagStrawHits(fhicl::ParameterSet const& pset) :
    _diag(pset.get<int>("diagLevel",0)),
    _debug(pset.get<int>("debugLevel",0)),
    _shTag(pset.get<string>("StrawHitCollection","makeSH")),
    _shpTag(pset.get<string>("StrawHitPositionCollection","MakeStereoHits")),
    _ccTag(pset.get<string>("CaloCluster","CaloClusterFromProtoCluster")),
    _usecc(pset.get<bool>("UseCaloCluster",false)),
    _useccxy(pset.get<bool>("UseCaloClusterXYPosition",false)),
    _ccmine(pset.get<float>("CaloClusterMinE",50.0)), // minimum energy required to a calorimeter cluster
    _ccminsize(pset.get<int>("CaloClusterMinSize")), // minimum calorimetercluster size
    _ccmindt(pset.get<float>("CaloClusterMinDt")), // minimum time residual of a stawhit from a calo cluster
    _ccmaxdt(pset.get<float>("CaloClusterMaxDt")), // maximum time residual of a stawhit from a calo cluster
    _pitchangle(pset.get<float>("pitchAngle")), //mean pitch angle of the particle we want to reconstruct
    _beta(pset.get<float>("beta")), //approximate relativistic beta of the particle we search
    _dtoffset(pset.get<float>("dtOffset")), //time offset of the strawhit form the calocluster
    _minE(pset.get<float>("minimumEnergy",0.0)), // Minimum deposited straw energy (MeV)
    _maxE(pset.get<float>("maximumEnergy",0.0035)), // MeV
    _ctE(pset.get<float>("crossTalkEnergy",0.007)), // MeV
    _ctMinT(pset.get<float>("crossTalkMinimumTime",-1)), // nsec
    _ctMaxT(pset.get<float>("crossTalkMaximumTime",100)), // nsec
    _minT(pset.get<float>("minimumTime",500)), // nsec
    _maxT(pset.get<float>("maximumTime",2000)), // nsec
    _minR(pset.get<float>("minimumRadius",395.0)), // mm
    _maxR(pset.get<float>("maximumRadius",650.0)), // mm
    _ttcalc(pset.get<fhicl::ParameterSet>("T0Calculator",fhicl::ParameterSet()))
  {
    produces<StrawHitFlagCollection>();
    _hmaxtdrift = 20.0; // ugly hack, replace with conditions access FIXME!
  }

  bool FlagStrawHits::findData(const art::Event& evt){
    _shcol  = 0;
    _shpcol = 0;
    _cccol  = 0;

    auto shH = evt.getValidHandle<StrawHitCollection>(_shTag);
    _shcol = shH.product();
    if(_shpTag.label() != ""){
      auto shpH = evt.getValidHandle<StrawHitPositionCollection>(_shpTag);
      _shpcol = shpH.product();
    }

    if(_usecc){
      auto ccH = evt.getValidHandle<CaloClusterCollection>(_ccTag);
      _cccol = ccH.product();
    }
    return _shcol != 0 && (_shpTag.label() == "" || _shpcol != 0) && (!_usecc || _cccol != 0);
  }

  void
  FlagStrawHits::produce(art::Event& event) {
    const Tracker& tracker = getTrackerOrThrow();
    const TTracker& tt = dynamic_cast<const TTracker&>(tracker);
    const Calorimeter* calorimeter; // cached pointer to the calorimeter geometry
    mu2e::GeomHandle<mu2e::Calorimeter> ch;
    calorimeter = ch.get();

    if ( _debug > 1 ) cout << "FlagStrawHits: produce() begin; event " << event.id().event() << endl;
    // find the data
    if(!findData(event)){
      throw cet::exception("RECO")<<"mu2e::FlagStrawHits: data missing or incomplete"<< endl;
      return;
    }

    size_t nsh = _shcol->size();
    // create the output collection
    unique_ptr<StrawHitFlagCollection> shfcol(new StrawHitFlagCollection);
    shfcol->reserve(nsh);
// A more efficient algorithm is to first find all big hits, then
// loop over only those in the float loop.  It avoids the search for indices FIXME!!

    size_t nplanes = tt.nPlanes();
    size_t npanels = tt.getPlane(0).nPanels();
    vector<vector<size_t> > hits_by_panel(nplanes*npanels,vector<size_t>());
   
   // select and sort hits by panel
    for (size_t ish=0;ish<nsh;++ish){
      StrawHit const& sh = _shcol->at(ish);
      const Straw& straw = tracker.getStraw( sh.strawIndex() );
      size_t iplane = straw.id().getPlane();
      size_t ipnl = straw.id().getPanel();
      size_t global_panel = ipnl + iplane*npanels;
      hits_by_panel[global_panel].push_back(ish);
    }

    vector<bool> is_ct_straw_neighbor(nsh,false);
    vector<bool> is_ct_straw_preamp(nsh,false);

    for (size_t ipanel=0;ipanel<nplanes*npanels;++ipanel){ // loop over panels
      // loop over hits in this panel
      for (size_t ish : hits_by_panel[ipanel]){
        StrawHit const& sh = _shcol->at(ish);
        const Straw& straw = tracker.getStraw( sh.strawIndex() );
        if (sh.energyDep() >= _ctE){
          for (size_t jsh : hits_by_panel[ipanel]){
            if (ish == jsh)
              continue;
            StrawHit const& sh2 = _shcol->at(jsh);
            if (sh2.time()-sh.time() > _ctMinT && sh2.time()-sh.time() < _ctMaxT){
              if (straw.isSamePreamp(sh2.strawIndex()))
                is_ct_straw_preamp[jsh] = true;
              if (straw.isNearestNeighbour(sh2.strawIndex()))
                is_ct_straw_neighbor[jsh] = true;
            }
          } // end loop over possible neighbors
        } // end if energy cut
      } // end loop over hits in panel
    } // end loop over panels

    // if we're using the calorimeter clusters find the time segments to search
    // vector<TimeRange> timeranges;
    // if(_usecc && _cccol != 0){
    //   timeranges.reserve(_cccol->size());
    //   for(auto const& cc : *_cccol ) {
    // 	if(goodCaloCluster(cc)){
    // 	// calculate the cluster time in the tracker reference
    // 	  float ctime = _ttcalc.caloClusterTime(cc);
    // 	  // add 1/2 the maximum drift time to center the selection
    // 	  ctime += _hmaxtdrift;
    // 	  // calculate the range around this that it's reasonable to search for tracks.  
    // 	  float tmin = ctime-_hctwin;
    // 	  float tmax = ctime+_hctwin;
    // 	  timeranges.push_back(make_pair(tmin,tmax));
    // 	}
    //   }
    // }
    // flags based on hit properties
    for(size_t ish=0;ish<nsh;++ish){
      StrawHit const& sh = _shcol->at(ish);
      StrawHitFlag flag;
      // merge with the position flag if that's present
      if(_shpcol != 0)flag.merge(_shpcol->at(ish).flag());
      if(sh.energyDep() > _minE && sh.energyDep() < _maxE)
        flag.merge(StrawHitFlag::energysel);
      if (is_ct_straw_neighbor[ish])
        flag.merge(StrawHitFlag::strawxtalk);
      if (is_ct_straw_preamp[ish])
        flag.merge(StrawHitFlag::elecxtalk);
      float time = std::min(sh.time(TrkTypes::cal),sh.time(TrkTypes::hv));
      if(time > _minT && time < _maxT)
        flag.merge(StrawHitFlag::timesel);
      // if(_usecc && goodTime(sh.time(),timeranges))
      // 	flag.merge(StrawHitFlag::calosel);
      if(_shpcol != 0){
        StrawHitPosition const& shp = _shpcol->at(ish);
        float rad = sqrtf(shp.pos().perp2());
        if(rad > _minR && rad < _maxR)
          flag.merge(StrawHitFlag::radsel);
      }
      shfcol->push_back(flag);
    }

    if(_usecc && _cccol != 0){
      int   ncl   =  _cccol->size();
      static const float pi(M_PI);
      static const float halfpi(pi/2.);
      static const float twopi(2*pi);
      float meanDriftTime = 1.25/0.06;// half straw tube radius / drift velocity; ugly hack, replace with conditions access FIXME!

      const CaloCluster*  cl(0);
      const StrawHit*     hit;
      const Straw*        straw;
      CLHEP::Hep3Vector   gpos, tpos;
      float              time, dt, tof, zstraw, cl_time;
      float              xcl, ycl, zcl, phicl;

      for (int ic=0; ic<ncl; ic++) {
	cl      = &_cccol->at(ic);
	
	if ( cl->energyDep() < _ccmine)        break;
	if ( (int(cl->size()) < _ccminsize) )  continue;
	
	cl_time = cl->time() + _dtoffset;
	//-----------------------------------------------------------------------------
	// convert cluster coordinates defined in the disk frame to the detector
	// coordinate system
	//-----------------------------------------------------------------------------
	gpos = calorimeter->geomUtil().diskToMu2e(cl->diskId(),cl->cog3Vector());
	tpos = calorimeter->geomUtil().mu2eToTracker(gpos);

	xcl     = tpos.x();
	ycl     = tpos.y();
	zcl     = tpos.z();
	phicl   = atan2f(ycl, xcl);

	//-----------------------------------------------------------------------------
	// record hits in time with each peak, and accept them if they have a minimum # of hits
	//-----------------------------------------------------------------------------

	for(size_t istr=0; istr<nsh;++istr) {
	  hit    = &_shcol->at(istr);
	  time   = hit->time();
	  straw  = &tt.getStraw(hit->strawIndex());
	  zstraw = straw->getMidPoint().z();
	  //-----------------------------------------------------------------------------
	  // estimate time-of-flight and calculate residual between the predicted and the hit times
	  // 2017-03-31 P.M.: this assumes electron (e^- or e^+), not muon
	  //-----------------------------------------------------------------------------
	  tof    = (zcl-zstraw)/sin(_pitchangle)/(CLHEP::c_light*_beta);
	  dt     = cl_time-(time+tof-meanDriftTime);
	  //-----------------------------------------------------------------------------
	  // 
	  //-----------------------------------------------------------------------------
	 

	  bool  caloCondition = (dt < _ccmaxdt) && (dt >= _ccmindt);

	  //use or not the calorimeter cluster XY coordinates to imrpove the selection?
	  if (_useccxy && (_shpcol != 0)){
	    float  phistraw = _shpcol->at(istr).pos().phi();
	    float  dphi     = phistraw - phicl; 
	    if (dphi >  pi) dphi -= twopi;
	    if (dphi < -pi) dphi += twopi;
	    caloCondition    = caloCondition && ( fabs(dphi) < halfpi );
	  }
	  
	  if ( caloCondition ) {
	    shfcol->at(istr).merge(StrawHitFlag::calopresel);
	  }
	  
	}//end loop over the straw hits
      }//end loop over the calorimeter clusters
    }

    event.put(move(shfcol));

  } // end FlagStrawHits::produce.

  bool FlagStrawHits::goodCaloCluster(CaloCluster const& cc) const {
    return cc.energyDep() > _ccmine;
  }

  // bool FlagStrawHits::goodTime(float hittime,vector<TimeRange> const& timeranges) const {
  //   bool retval(false);
  //   for(auto const& trange : timeranges ){
  //     if(hittime > trange.first && hittime < trange.second){
  // 	retval = true;
  // 	break;
  //     }
  //   }
  //   return retval;
  // }
 
} // end namespace mu2e

using mu2e::FlagStrawHits;
DEFINE_ART_MODULE(FlagStrawHits)

