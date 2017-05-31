//
// A module to create simple stereo hits out of StrawHits.  This can work
// with either tracker.  StrawHit selection is done by flagging in an upstream module
//
// $Id: MakeStereoHits_module.cc,v 1.23 2014/09/18 08:42:47 brownd Exp $
// $Author: brownd $
// $Date: 2014/09/18 08:42:47 $
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
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/StereoHit.hh"
#include "Mu2eUtilities/inc/FMVATools.hh"

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
#include <float.h>
#include "CalPatRec/inc/THackData.hh"

//#include "PerfLib/inc/perflib.hh"
//perf::PerfStats g_perf("MakeStereoHits 100") ;

using namespace std;

#include "PerfLib/inc/perflib.hh"
using perf::PerfStats;
using perf::metrics_t;

namespace mu2e {
  // struct to hold MVA input  
  // struct to hold MVA input
  struct StereoMVA {
    std::vector <Double_t> _pars;
    Double_t& _dt; // time diff between hits
    Double_t& _chisq; // chisq of time division information 
    Double_t& _rho;  // transverse radius of stereo position
    Double_t& _ndof; // number of degrees of freedom of time division chisquared: either 0, 1, or 2
    StereoMVA() : _pars(4,0.0),_dt(_pars[0]),_chisq(_pars[1]),_rho(_pars[2]),_ndof(_pars[3]){}
  };

  class MakeStereoHits : public art::EDProducer {

    typedef vector<size_t> PanelHits;
    typedef vector<PanelHits> StationPanels;
    typedef vector< StationPanels > TrackerStations;

  public:
    explicit MakeStereoHits(fhicl::ParameterSet const& pset);

    void produce( art::Event& e);
    virtual void beginJob();
    virtual void beginRun(art::Run & run);

  private:

    // configuration
    int _debug;
    // event object Tags
    art::InputTag   _shTag;
    art::InputTag   _shpTag;
    art::InputTag   _shfTag;
    // input event collection cache
    const StrawHitCollection*             _shcol;
    const StrawHitPositionCollection*     _shpcol;
    const StrawHitFlagCollection*	  _shfcol;
    // Parameters
    StrawHitFlag _shsel; // flag selection
    StrawHitFlag _shmask; // flag anti-selection 
    double _maxDt; // maximum time separation between hits
    double _maxDE; // maximum deposited energy deference: this excludes inconsistent hits
    double _maxDZ; // maximum longitudinal separation
    double _maxDPerp; // maximum transverse separation
    double _minDdot; // minimum dot product of straw directions
    double _minDL; // minimum distance from end of active straw;
    double _maxChi; // maximum # of TimeDivision sigmas past the active edge of a wire to allow making stereo hits
    double _maxChisq; // maximum # of TimeDivision consistency chisquared to allow making stereo hits
    double _minMVA; // minimum MVA output
    double _wres; // resolution to assign along the wire
    bool _bestpair; // require both hits agree as to the best pair
    bool _writepairs; // write out the stereo pairs
    // MVA
    MVATools _mvatool;
    StereoMVA _vmva; // input variables to TMVA for stereo selection
    // for optimized Stereo Hit finding
    size_t _nsta; // number of stations
    size_t _npnl; // number of panels in a panel
    vector <vector<size_t> >_pover;            // list of panels overlapping a given panel
    // helper functions
    void genMap();    // function to generate panel list from tracker geometry
    double longRes(StereoHit const& sthit) const; // longitudinal resolution
    bool betterPair(StereoHit const& newpair, StereoHit const& oldpair) const;
    bool findData(art::Event& event);
    void reportInsertedHits(TrackerStations const & stax);
 };

  MakeStereoHits::MakeStereoHits(fhicl::ParameterSet const& pset) :
    // Parameters
    _debug(pset.get<int>("debugLevel",0)),
    _shTag		(pset.get<art::InputTag>("StrawHitCollection","makeSH")),
    _shpTag		(pset.get<art::InputTag>("StrawHitPositionCollection","MakeStrawHitPositions")),
    _shfTag	       (pset.get<art::InputTag>("StrawHitFlagCollection","FSHPreStereo")),
    _shsel(pset.get<vector<string> >("StrawHitSelectionBits",vector<string>{"EnergySelection","TimeSelection"} )),
    _shmask(pset.get<vector<string> >("StrawHitMaskBits",vector<string>{} )),
    _maxDt(pset.get<double>("maxDt",40.0)), // nsec
    _maxDE(pset.get<double>("maxDE",0.99)), // dimensionless, from 0 to 1
    _maxDZ(pset.get<double>("maxDZ",1000.)), // mm, maximum longitudinal distance between straws
    _maxDPerp(pset.get<double>("maxDPerp",500.)), // mm, maximum perpendicular distance between time-division points
    _minDdot(pset.get<double>("minDdot",0.6)), // minimum angle between straws
    _minDL(pset.get<double>("minDL",-20.0)), // extent along straw
    _maxChisq(pset.get<double>("maxChisquared",80.0)), // position matching
    _minMVA(pset.get<double>("minMVA",0.6)), // MVA cut
    _wres(pset.get<double>("LongitudinalResolution",20.0)), // estimated resolution of stereo reco
    _bestpair(pset.get<bool>("BestStereoPair",false)),
    _writepairs(pset.get<bool>("WriteStereoPairs",true)),
    _mvatool(pset.get<fhicl::ParameterSet>("MVATool",fhicl::ParameterSet()))
  {
    _maxChi = sqrt(_maxChisq);
    // Tell the framework what we make.
    if(_writepairs)produces<StereoHitCollection>();
    produces<StrawHitPositionCollection>();
  }

  void MakeStereoHits::beginJob(){
  // initialize MVA
    _mvatool.initMVA();
    if(_debug > 0){
      cout << "MakeStereoHits MVA parameters: " << endl;
      _mvatool.showMVA();
    }
  }

  void MakeStereoHits::beginRun(art::Run & run){
    genMap();
  }

void MakeStereoHits::reportInsertedHits(TrackerStations const & stax){
      if(!stax.empty())	{
	if(!stax[0].empty())
	  cout << "Built PanelHits size = " << stax[0][0].size() << endl;
	cout << "Built StationPanels size = " << stax[0].size() << endl;
      }
      cout << "Built TrackerStations size = " << stax.size() << endl;
      for( size_t station =0; station < stax.size(); ++station) {                      // loop over stations
	for(size_t ipnl = 0; ipnl < stax[station].size(); ++ipnl) {   // loop over panels in this station
	  for(size_t ish : stax[station][ipnl]){ 	   
	    cout << "Inserted hit " << ish << " into station " << station << " panel " << ipnl << endl;
	  }
	}
      }
}

struct free_delete{
    void operator()(void* x) { free(x); }
};

  void MakeStereoHits::produce(art::Event& event) {
  //   g_perf.read_begin_counters_inlined();

    // Get a reference to T trackers
    const Tracker& tracker = getTrackerOrThrow();
    const TTracker& tt = dynamic_cast<const TTracker&>(tracker);

    // find event data
    if( !findData(event) ){
      throw cet::exception("RECO")<<"mu2e::MakeStereoHits: collection missing or sizes don't match " << endl; 
    }

    // create the stereo hits
    StereoHitCollection stereohits;
    size_t nsh = _shcol->size();
    stereohits.reserve(3*nsh);
    // index to best stereo hit for a given straw hit
    vector<int> ibest(nsh,-1);
    // select and sort hits by panel
    size_t nres = max(size_t(32),nsh/10);

    PanelHits hpstax;
    hpstax.reserve(nres);
    StationPanels pstax(2*_npnl,hpstax);
    TrackerStations stax(_nsta,pstax);

    for(size_t ish=0;ish<nsh;++ish){
      // merge the flags
      StrawHitFlag shf = _shfcol->at(ish);
      shf.merge(_shpcol->at(ish).flag());
    // select hits based on flag
      if(shf.hasAllProperties(_shsel) && (!shf.hasAnyProperty(_shmask)) ){
	StrawHit const& hit = _shcol->at(ish);
	Straw const& straw = tt.getStraw(hit.strawIndex());
        size_t iplane = straw.id().getPlane();
        size_t ipnl = straw.id().getPanel();
        size_t station = iplane/2;
	// define a 'global' panel for the station.  This changes sign with odd-even stations
	size_t jpnl;
	if(station%2==0)
	  jpnl = ipnl + (iplane%2)*_npnl;
	else
	  jpnl = ipnl + (1-iplane%2)*_npnl;
	        
        stax[station][jpnl].push_back(ish);
      }
    }
    
    if(_debug > 2) 
      reportInsertedHits(stax);

    // prepair hitpairs    
    size_t nhitpairs=0;
        
    for(StationPanels const& pstax : stax){ // loop over stations
      auto pstax_size = pstax.size();      
      for(size_t ipnl = 0; ipnl < pstax_size; ++ipnl) {             // loop over panels in this station	
	if( _pover[ipnl].empty() || pstax[ipnl].empty() ) continue;	
	  for( size_t jpnl : _pover[ipnl]){                        // loop over overlapping panels	      
	      nhitpairs+=pstax[ipnl].size()*pstax[jpnl].size();
	  }
	}
      }
    
    std::unique_ptr<size_t, free_delete>buffer((size_t*)malloc(2*nhitpairs*sizeof(size_t)));
    size_t*const shprcol=buffer.get();
        
    size_t idx=0;

    for(StationPanels const& pstax : stax){  // loop over stations
      auto pstax_size = pstax.size();
      for(size_t ipnl = 0; ipnl < pstax_size; ++ipnl) {             // loop over panels in this station	
	if( (!_pover[ipnl].empty()) && (!pstax[ipnl].empty()) ){
	  for( size_t jpnl : _pover[ipnl]){                        // loop over overlapping panels
	    if( !pstax[jpnl].empty() ){
	      for(size_t ish : pstax[ipnl]){                                    // loop over selected hits in the first panel
		for(size_t jsh : pstax[jpnl]){
		    shprcol[idx] =ish;
		    shprcol[idx+1] =jsh;
		    idx+=2;
		}
	      }	      
	    }
	  }
	}
      }
    }

   // loop over hitpairs
  for(size_t ishpr=0;ishpr < idx;ishpr+=2) {
 		  auto ish=shprcol[ishpr];
		  auto jsh=shprcol[ishpr+1];
		  
		  StrawHit const& sh1 = _shcol->at(ish);
		  StrawHit const& sh2 = _shcol->at(jsh);

   		  double dt = fabs(sh1.time()-sh2.time()); //acceptance  0.11
		  if(dt >= _maxDt )// hits are close in time
		    continue;
		  Straw const& straw1 = tt.getStraw(sh1.strawIndex());
		  Straw const& straw2 = tt.getStraw(sh2.strawIndex());
		  double ddot = straw1.direction().dot(straw2.direction()); //acceptance 0.35
		  if( ddot <= _minDdot )// hits are close in time
		    continue;
		  StrawHitPosition const& shp1 = _shpcol->at(ish);
		  StrawHitPosition const& shp2 = _shpcol->at(jsh);		  

		  CLHEP::Hep3Vector dp = shp1.pos()-shp2.pos();
		  double dperp = dp.perp(); //acceptance 0.47

		  if (dperp >= _maxDPerp)
		    continue;
		  double de = min((float)1.0, fabs((sh1.energyDep() - sh2.energyDep())/(sh1.energyDep()+sh2.energyDep())));
		  if (de >= _maxDE)
		    continue;
		  //double dz = fabs(dp.z()); //acceptance 1.0
		  //PanelId::isep sep = straw1.id().getPanelId().separation(straw2.id().getPanelId()); //acceptance 1.0
		  if( //sep != PanelId::same && sep < PanelId::apart // hits are in the same station but not the same panel
		      ddot > _minDdot // negative crosings are in opposite quadrants
		      && dt < _maxDt // hits are close in time
		      && de < _maxDE   // deposited energy is roughly consistent (should compare dE/dx but have no path yet!)
		     // && dz < _maxDZ // longitudinal separation isn't too big
		      && dperp < _maxDPerp
		       ) { // transverse separation isn't too big
		    // tentative stereo hit: this solves for the POCA
		    StereoHit sth(*_shcol,tt,ish,jsh);
		    double dl1 = straw1.getDetail().activeHalfLength()-fabs(sth.wdist1());
		    double dl2 = straw2.getDetail().activeHalfLength()-fabs(sth.wdist2());
		    if( dl1 > _minDL && dl2 > _minDL) {
		      // stereo point is inside the active length
		      // compute difference between stereo points and TD prediction
		      double chi1(0.0), chi2(0.0);
		      unsigned ndof(0);

		      bool tdiv1 = shp1.flag().hasAllProperties(StrawHitFlag::tdiv);
		      if(tdiv1){
			chi1 = (shp1.wireDist()-sth.wdist1())/shp1.posRes(StrawHitPosition::wire);
			++ndof;
		      }

		      bool tdiv2 = shp2.flag().hasAllProperties(StrawHitFlag::tdiv);
		      if(tdiv2){
			chi2 = (shp2.wireDist()-sth.wdist2())/shp2.posRes(StrawHitPosition::wire);
			++ndof;
		      }
		      if(fabs(chi1) <_maxChi && fabs(chi2) < _maxChi) {
			// compute chisquared
			double chisq = chi1*chi1+chi2*chi2; 
			if(chisq < _maxChisq){
			  sth.setChisquared(chisq);
			  // compute MVA
			  _vmva._dt = dt;
			  _vmva._chisq = chisq;
			  _vmva._rho = sth.pos().perp();
			  _vmva._ndof = ndof;
			  double mvaout = _mvatool.evalMVA(_vmva._pars);
			  if(mvaout > _minMVA){
			  // valid stereo pairing.  Record it
			    sth.setMVAOut(mvaout);
			    stereohits.push_back(sth);
			    size_t isth = stereohits.size()-1;
			   //  See if this is better than the existing pairs 
			    if(ibest[ish] < 0 || betterPair(sth,stereohits.at(ibest[ish])))
			      ibest[ish] = isth;
			    if(ibest[jsh] < 0 || betterPair(sth,stereohits.at(ibest[jsh])))
			      ibest[jsh] = isth;
			  } // good MVA
			} // good combined chisquared
		      } // good individual chis
		    } // good longitudinal separation
		  } // good initial hit compatibility
}
    // now, overwrite the positions for those hits which have stereo
    // copy the input StrawHitPosition collection
    unique_ptr<StrawHitPositionCollection> shpcol(new StrawHitPositionCollection(*_shpcol));
    for(size_t ish=0; ish<nsh;++ish){
      if(ibest[ish] > 0){
	StereoHit const& sthit = stereohits.at(ibest[ish]);
	// optionally require this be the best pairing for both hits
	if( (!_bestpair) || ibest[sthit.hitIndex1()] == ibest[sthit.hitIndex2()]){
	  StrawHitPosition& shpos = shpcol->at(ish);
	  shpos._flag.merge(StrawHitFlag::stereo);
	  shpos._stindex = ibest[ish];
	  shpos._pos = sthit.pos();
	  shpos._wres = longRes(sthit);
	}
      }
    }
    // if requested, put the stereo hits themselves into the event
    if(_writepairs){
      unique_ptr<StereoHitCollection> sthits(new StereoHitCollection(move(stereohits)));
      event.put(move(sthits));
    }
    event.put(move(shpcol));
   // g_perf.read_end_counters_inlined();
  } // end MakeStereoHits::produce.

// estimate the resolution on the stereo hit position projection along a wire direction
  double MakeStereoHits::longRes(StereoHit const& sthit) const {
    return _wres; // this should be a real calculation FIXME!!!
  }

  bool MakeStereoHits::betterPair(StereoHit const& newpair, StereoHit const& oldpair) const {
    // choose the best pair as:
    // 1) take the pair with the minimum plane separation
    // 2) otherwise, take the pair with the maximum MVA output
    if(newpair.panelSeparation() == oldpair.panelSeparation() ) {
      return newpair.mvaout() > oldpair.mvaout();
    } else {
      return newpair.panelSeparation() < oldpair.panelSeparation();
    }
  }

  void MakeStereoHits::genMap() {
    // setup panel map
    const Tracker& tracker = getTrackerOrThrow();
    const TTracker& tt = dynamic_cast<const TTracker&>(tracker);
    // establish sizes
    _nsta = tt.nPlanes()/2;
    _npnl = tt.getPlane(0).nPanels(); // number of panels in a plane
    // find the phi extent of the longest straw
    Straw const& straw = tt.getStraw(StrawId(0,0,0,0));
    double phi0 = (straw.getMidPoint()-straw.getHalfLength()*straw.getDirection()).phi();
    double phi1 = (straw.getMidPoint()+straw.getHalfLength()*straw.getDirection()).phi();
    double lophi = min(phi0,phi1);
    double hiphi = max(phi0,phi1);
    double phiwidth = hiphi-lophi;
    if(phiwidth>M_PI)phiwidth =2*M_PI-phiwidth;
    // loop over stations and see whether the phi ranges of the straws overlap, giving
    // a possibility for stereo hits
    std::vector<double> panphi(2*_npnl);

    for(int ipla=0;ipla<2;++ipla){
      Plane const& plane = tt.getPlane(ipla);
      for(int ipan=0;ipan<plane.nPanels();++ipan){
	Panel const& panel = plane.getPanel(ipan);
	// expand to station-wide 'panel' number.  This changes sign with station
	int jpan = ipan + (ipla%2)*_npnl;
	panphi[jpan] = panel.straw0MidPoint().phi();
      }
    }

    if(_debug>0){
      cout << "panel phi width = " << phiwidth;
      cout << "panel phi positions = ";
      for(auto pphi : panphi)
	cout << pphi << " ";
      cout << endl;
    }

    _pover = vector<vector<size_t>>(2*_npnl);
    if(_debug< 10){
      for(size_t iphi = 0;iphi<2*_npnl;++iphi){
	double phi = panphi[iphi];
	for(size_t jphi=iphi+1;jphi<2*_npnl;++jphi){
	  double dphi = fabs(phi - panphi[jphi]);
	  if(dphi > M_PI)dphi = 2*M_PI-dphi;
	  if(dphi < phiwidth)_pover[iphi].push_back(jphi);
	}
      }
    } else {
      // for testing, assume every panel overlaps
      for(size_t iphi = 0;iphi<2*_npnl;++iphi){
	for(size_t jphi=iphi+1;jphi<2*_npnl;++jphi){
	  _pover[iphi].push_back(jphi);
	}
      }
    }

    if(_debug >0){
      for(size_t ipnl=0;ipnl<2*_npnl;++ipnl){
	cout << "Panel " << ipnl << " Overlaps with panel ";
	for(auto jpnl : _pover[ipnl]){
	  cout << jpnl << " ";
	}
	cout << endl;
      }
    }
  }

  bool MakeStereoHits::findData(art::Event& evt){
    _shcol = 0; _shpcol = 0;  _shfcol = 0;
    auto shH = evt.getValidHandle<StrawHitCollection>(_shTag);
    _shcol = shH.product();
    auto shpH = evt.getValidHandle<StrawHitPositionCollection>(_shpTag);
    _shpcol = shpH.product();
    auto shfH = evt.getValidHandle<StrawHitFlagCollection>(_shfTag);
    _shfcol = shfH.product();
    return _shcol != 0 && _shpcol != 0 && _shpcol->size() == _shcol->size();
  }


} // end namespace mu2e

using mu2e::MakeStereoHits;
DEFINE_ART_MODULE(MakeStereoHits)

