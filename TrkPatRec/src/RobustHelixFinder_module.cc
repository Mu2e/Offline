//
// TTracker Pattern Recognition based on Robust Helix Fit
//
// $Id: RobustHelixFinder_module.cc,v 1.2 2014/08/30 12:19:38 tassiell Exp $
// $Author: tassiell $
// $Date: 2014/08/30 12:19:38 $
//
// Original author D. Brown and G. Tassielli
//

// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
// utilities
#include "GeneralUtilities/inc/Angles.hh"
#include "Mu2eUtilities/inc/MVATools.hh"
/// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/TimeClusterCollection.hh"
#include "RecoDataProducts/inc/HelixSeedCollection.hh"
#include "RecoDataProducts/inc/TrkFitFlag.hh"
// tracking
#include "TrkReco/inc/TrkTimeCalculator.hh"
// BaBar
#include "BTrk/BaBar/BaBar.hh"
#include "TrkReco/inc/TrkDef.hh"
#include "BTrkData/inc/TrkStrawHit.hh"
#include "TrkReco/inc/RobustHelixFit.hh"
//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"
// boost
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/accumulators/statistics/median.hpp>
// root
#include "TH1F.h"
// C++
#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <utility>
#include <functional>
#include <float.h>
#include <vector>
#include <set>
#include <map>
using namespace std; 
using CLHEP::Hep3Vector;
using CLHEP::HepVector;
using CLHEP::HepSymMatrix;
using namespace boost::accumulators;

namespace mu2e 
{

  // struct to hold MVA input
  struct HelixHitMVA {
    std::vector <Double_t> _pars;
    Double_t& _dtrans; // distance from hit to helix perp to the wrire
    Double_t& _dwire; // distance from hit to helix along the wrire
    Double_t& _chisq; // chisq of spatial information, using average errors
    Double_t& _dt;  // time difference of hit WRT average
    Double_t& _drho; // hit transverse radius minus helix radius
    Double_t& _dphi; // hit azimuth minus helix azimuth (at the hit z)
    Double_t& _rwdot; // dot product between circle radial direction and wire direction
    Double_t& _hrho; // helix transverse radius (at the hit z)
    Double_t& _hhrho; // hit transverse radius
    HelixHitMVA() : _pars(9,0.0),_dtrans(_pars[0]),_dwire(_pars[1]),_chisq(_pars[2]),_dt(_pars[3]),
     _drho(_pars[4]),_dphi(_pars[5]),_rwdot(_pars[6]),_hrho(_pars[7]),_hhrho(_pars[8]){}
  };

  class RobustHelixFinder : public art::EDProducer
  {
  public:
    explicit RobustHelixFinder(fhicl::ParameterSet const&);
    virtual ~RobustHelixFinder();
    virtual void beginJob();
    virtual void produce(art::Event& event ); 
  private:
    unsigned                           _iev;

    // configuration parameters
    int                                 _diag,_debug;
    int                                 _printfreq;
    bool				_prefilter; // prefilter hits based on sector
    unsigned				_minnhit; // minimum # of hits to work with
    double				_maxdr; // maximum hit-helix radius difference
    double				_maxrpull; // maximum hit-helix radius difference pull
    double				_maxphisep; // maximum separation in global azimuth of hits
    TrkFitFlag				_saveflag; // write out all helices that satisfy these flags
    unsigned				_maxniter;  // maximum # of iterations over outlier filtering + fitting
    double				_cradres; // average center resolution along center position (mm)
    double				_cperpres; // average center resolution perp to center position (mm)
    double				_radres; // average radial resolution for circle fit (mm)
    double				_phires; // average azimuthal resolution on circle (rad)
    double				_maxdwire; // outlier cut on distance between hit and helix along wire
    double				_maxdtrans; // outlier cut on distance between hit and helix perp to wire
    double				_maxchisq; // outlier cut on chisquared
    double				_maxrwdot[2]; // outlier cut on angle between radial direction and wire: smaller is better
    
    bool				_usemva; // use MVA to cut outliers
    double				_minmva; // outlier cut on MVA

    // input object tags
    art::InputTag			_shTag;
    art::InputTag			_shpTag;
    art::InputTag			_shfTag;
    art::InputTag			_tcTag;
    // output label
    std::string                        _trackseed;

    // initial hit selection
    StrawHitFlag  _hsel, _hbkg;

  // MVA for iteratively cleaning hits from the helix
    MVATools _mvatool;
    HelixHitMVA _vmva; // input variables to TMVA for filtering hits
    // Histograms
    TH1F* _niter, *_nitermva;

    // cache of event objects
    const StrawHitCollection*  _shcol;
    const StrawHitPositionCollection*  _shpcol;
    const StrawHitFlagCollection*      _shfcol;
    const TimeClusterCollection*       _tccol;

    // robust helix fitter
    RobustHelixFit                     _hfit;
    TrkTimeCalculator			_ttcalc;
    double _t0shift;   // adhoc t0 shift
    bool			      _parcor; // correct the parameters
    std::vector<double> _crcorr, _rcorr; // corrections for center radius and radius
    // helper functions
    bool findData           (const art::Event& e);
    void prefilterHits(HelixSeed& hseed); // rough filtering based on hits being in 1/2 the detector
    bool filterCircleHits(HelixSeed& hseed); // return value tells if any hits changed state
    bool filterHits(HelixSeed& hseed); // return value tells if any hits changed state
    void fillMVA(HelixSeed& hseed); // return value tells if any hits changed state
    bool filterHitsMVA(HelixSeed& hseed); // return value tells if any hits changed state
    void updateT0(HelixSeed& hseed); // update T0 value based on current good hits
    void correctParameters(RobustHelix& helix);// correct the helix parameters to MC truth using linear relations and correlations
    unsigned hitCount(HelixSeed const& hseed);
    void fitHelix(HelixSeed& hseed);
 
  };

  RobustHelixFinder::RobustHelixFinder(fhicl::ParameterSet const& pset) :
    _diag        (pset.get<int>("diagLevel",0)),
    _debug       (pset.get<int>("debugLevel",0)),
    _printfreq   (pset.get<int>("printFrequency",101)),
    _prefilter   (pset.get<bool>("PrefilterHits",true)),
    _minnhit	 (pset.get<unsigned>("minNHit",5)),
    _maxdr	 (pset.get<double>("MaxRadiusDiff",40.0)), // mm
    _maxrpull	 (pset.get<double>("MaxRPull",10.0)), // unitless
    _maxphisep	 (pset.get<double>("MaxPhiHitSeparation",1.0)),
    _saveflag    (pset.get<vector<string> >("SaveHelixFlag",vector<string>{"HelixOK"})),
    _maxniter    (pset.get<unsigned>("MaxIterations",10)), // iterations over outlier removal
    _cradres	 (pset.get<double>("CenterRadialResolution",12.0)),
    _cperpres	 (pset.get<double>("CenterPerpResolution",12.0)),
    _radres	 (pset.get<double>("RadiusResolution",10.0)),
    _phires	 (pset.get<double>("AzimuthREsolution",0.1)),
    _maxdwire    (pset.get<double>("MaxWireDistance",200.0)), // max distance along wire
    _maxdtrans   (pset.get<double>("MaxTransDistance",100.0)), // max distance perp to wire (and z)
    _maxchisq    (pset.get<double>("MaxChisquared",100.0)), // max chisquared
    _usemva      (pset.get<bool>("UseHitMVA",false)),
    _minmva      (pset.get<double> ("MinMVA",0.2)), // min MVA output to define an outlier
    _shTag	 (pset.get<art::InputTag>("StrawHitCollection","makeSH")),
    _shpTag	 (pset.get<art::InputTag>("StrawHitPositionCollection","MakeStereoHits")),
    _shfTag	 (pset.get<art::InputTag>("StrawHitFlagCollection","TimeClusterFinder")),
    _tcTag	 (pset.get<art::InputTag>("TimeClusterCollection","TimeClusterFinder")),
    _trackseed   (pset.get<string>("HelixSeedCollectionLabel","TimeClusterFinder")),
    _hsel        (pset.get<std::vector<std::string> >("HitSelectionBits")),
    _hbkg        (pset.get<vector<string> >("HitBackgroundBits",vector<string>{"DeltaRay","Isolated"})),
    _mvatool     (pset.get<fhicl::ParameterSet>("HelixHitMVA",fhicl::ParameterSet())),
    _hfit        (pset.get<fhicl::ParameterSet>("RobustHelixFit",fhicl::ParameterSet())),
    _ttcalc      (pset.get<fhicl::ParameterSet>("T0Calculator",fhicl::ParameterSet())),
    _t0shift     (pset.get<double>("T0Shift",0.0)),
    _parcor      (pset.get<bool>("CorrectHelixParameters",false)),
    _crcorr      (pset.get<vector<double> >("CenterRadiusCorrection")),
    _rcorr       (pset.get<vector<double> >("RadiusCorrection"))
  {
    _maxrwdot[0] = pset.get<double>("MaxStereoRWDot",1.0);
    _maxrwdot[1] = pset.get<double>("MaxNonStereoRWDot",0.9);
    produces<HelixSeedCollection>();
  }

  RobustHelixFinder::~RobustHelixFinder(){}

  void RobustHelixFinder::beginJob() {
  // initialize MVA
    _mvatool.initMVA();
    if(_debug > 0){
      cout << "RobustHeilxFinder MVA parameters: " << endl;
      _mvatool.showMVA();
    }
    if(_diag > 0){
      art::ServiceHandle<art::TFileService> tfs;
      _niter = tfs->make<TH1F>( "nitever" , "Number of Fit Iteraions",201,-0.5,200.5);
      _nitermva = tfs->make<TH1F>( "nitermva" , "Number of MVA Fit Iteraions",201,-0.5,200.5);
    }
  }

  void RobustHelixFinder::produce(art::Event& event ) {

    // create output
    unique_ptr<HelixSeedCollection> outseeds(new HelixSeedCollection);
    // event printout
    _iev=event.id().event();
    // find the data
    if(!findData(event)){
      throw cet::exception("RECO")<<"mu2e::RobustHelixFinder: data missing or incomplete"<< endl;
    }

    for(auto itclust=_tccol->begin(); itclust != _tccol->end(); ++ itclust) {
      auto const& tclust = *itclust;
    // build an empty HelixSeed 
      HelixSeed hseed;
      // copy in the t0 and cluster
      hseed._t0 = tclust._t0;
      hseed._caloCluster = tclust._caloCluster;
      // create a ptr to the time cluster
      size_t index = std::distance(_tccol->begin(),itclust);
      auto tcH = event.getValidHandle<TimeClusterCollection>(_tcTag);
      hseed._timeCluster = art::Ptr<TimeCluster>(tcH,index);
      unsigned niter(0);
      // loop over hits in this time cluster and select  hits with good 3-d position information
      std::vector<StrawHitIndex> goodhits;
      for(auto const& ind : tclust._strawHitIdxs) {
	if(_shfcol->at(ind).hasAnyProperty(_hsel) && !_shfcol->at(ind).hasAnyProperty(_hbkg))
	  goodhits.push_back(ind);
      }
      // create helix seed hits from the straw hit positions 
      for(auto idx : goodhits ) {
	HelixHit hhit(_shpcol->at(idx),idx);
	hhit._flag.clear(StrawHitFlag::resolvedphi);
	// merge in the other flags
	hhit._flag.merge(_shfcol->at(idx));
	hseed._hhits.push_back(hhit);
      }
      // prefilter hits
      if(_prefilter)prefilterHits(hseed);
      // check we have enough hits
      if(hitCount(hseed) >= _minnhit){
	hseed._status.merge(TrkFitFlag::hitsOK);
	// fit the helix with simple hit filtering
      	fitHelix(hseed);
	if(_usemva) {
      // repeat, using MVA filtering
	  niter = 0;
	  bool changed = true;
	  while(hseed._status.hasAllProperties(TrkFitFlag::helixOK)  && niter < _maxniter && changed) {
	    // fill hit MVA information
	    fillMVA(hseed);
	    // MVA filtering is sensitive to both radial and phi-Z
	    changed = filterHitsMVA(hseed);
	    _hfit.fitHelix(hseed);
	    // update t0 as that's part of the MVA
	    if(changed)updateT0(hseed);
	    ++niter;
	  }
	  if(_diag > 0)_nitermva->Fill(niter);
	}
      } else 
      // simply fill the MVA information for diagnostics
	fillMVA(hseed);
      // test fit status
      if(hseed.status().hasAllProperties(_saveflag)){
	if(niter < _maxniter)hseed._status.merge(TrkFitFlag::helixConverged);
// correct the parameters if requested
	if(_parcor)correctParameters(hseed._helix);
// update t0
	updateT0(hseed);
// save this helix
	outseeds->push_back(hseed);
	if(_debug > 1) cout << "Found helix with fit \n" << hseed._helix << endl;
      } else if (_debug > 1) cout << "Found helix without fit \n" << hseed._helix << endl;
    }

    if (_debug>0 && (_iev%_printfreq)==0) {
      cout<<"event "<<_iev<<" tot N hit "<<_shfcol->size()<<" N tracks seed found "<<outseeds->size()
	       <<" N time peaks "<<_tccol->size()<<endl;
    }

    event.put(std::move(outseeds));
  }

  void RobustHelixFinder::fillMVA(HelixSeed& hseed) {
    RobustHelix& helix = hseed._helix;
    HelixHitCollection& hhits = hseed._hhits;
    static Hep3Vector zaxis(0.0,0.0,1.0); // unit in z direction
// loop over hits
    for(auto& hhit : hhits) {
      // directions
      Hep3Vector const& wdir = hhit.wdir();
      Hep3Vector wtdir = zaxis.cross(wdir); // transverse direction to the wire
      Hep3Vector cvec = (hhit.pos() - helix.center()).perpPart(); // direction from the circle center to the hit
      Hep3Vector cdir = cvec.unit(); // direction from the circle center to the hit
      Hep3Vector cperp = zaxis.cross(cdir); // direction perp to the radius
      // positions
      Hep3Vector hpos = hhit.pos(); // this sets the z position to the hit z
      helix.position(hpos); // this computes the helix expectation at that z
      Hep3Vector dh = hhit.pos() - hpos; // this is the vector between them
      // fill MVA struct
      _vmva._dtrans = fabs(dh.dot(wtdir)); // transverse projection
      _vmva._dwire = fabs(dh.dot(wdir)); // projection along wire direction
      _vmva._drho = fabs(cvec.mag() - helix.radius()); // radius difference
      _vmva._dphi = fabs(hhit._phi - helix.circleAzimuth(hhit.pos().z())); // azimuth difference WRT circle center
      _vmva._hhrho = cvec.mag(); // hit transverse radius WRT circle center
      _vmva._hrho = hpos.perp();  // hit detector transverse radius
      _vmva._rwdot = fabs(wdir.dot(cdir)); // compare directions of radius and wire
      // compute the total resolution including hit and helix parameters first along the wire
      double wres2 = std::pow(hhit.posRes(StrawHitPosition::wire),(int)2) +
	std::pow(_cradres*cdir.dot(wdir),(int)2) +
	std::pow(_cperpres*cperp.dot(wdir),(int)2);
      // transverse to the wires
      double wtres2 = std::pow(hhit.posRes(StrawHitPosition::trans),(int)2) +
	std::pow(_cradres*cdir.dot(wtdir),(int)2) +
	std::pow(_cperpres*cperp.dot(wtdir),(int)2);
      // create a chisquared from these
      _vmva._chisq = sqrt( _vmva._dwire*_vmva._dwire/wres2 + _vmva._dtrans*_vmva._dtrans/wtres2 );
      // time difference
      _vmva._dt = _shcol->at(hhit._shidx).time() - hseed._t0.t0();
      // compute the MVA and store it
      hhit._hqual = _mvatool.evalMVA(_vmva._pars);
    }
  }
  
  bool RobustHelixFinder::filterHitsMVA(HelixSeed& hseed) {
    HelixHitCollection& hhits = hseed._hhits;
    bool changed(false);
    static StrawHitFlag outlier(StrawHitFlag::outlier);
// loop over hits
    for(auto& hhit : hhits) {
      bool oldout = !hhit._flag.hasAnyProperty(outlier);
      if(hhit._hqual < _minmva ) {
	// outlier hit flag it
	hhit._flag.merge(outlier);
      } else {
	// clear the hit in case it was formerly an outlier
	hhit._flag.clear(outlier);
      }
      changed |= oldout != hhit._flag.hasAnyProperty(outlier);
    }
    return changed;
  }

  bool RobustHelixFinder::filterCircleHits(HelixSeed& hseed) {
    RobustHelix& helix = hseed._helix;
    HelixHitCollection& hhits = hseed._hhits;
    bool changed(false);
    static StrawHitFlag outlier(StrawHitFlag::outlier);
    static Hep3Vector zaxis(0.0,0.0,1.0); // unit in z direction
// loop over hits
    for(auto& hhit : hhits) {
      bool oldout = !hhit._flag.hasAnyProperty(outlier);
      // directions
      Hep3Vector const& wdir = hhit.wdir();
      Hep3Vector cvec = (hhit.pos() - helix.center()).perpPart(); // direction from the circle center to the hit
      Hep3Vector cdir = cvec.unit(); // direction from the circle center to the hit
      double rwdot = wdir.dot(cdir); // compare directions of radius and wire
      // compute radial difference and pull
      double dr = cvec.mag()-helix.radius();
      double werr = hhit.posRes(StrawHitPosition::wire);
      // the resolution is dominated
      double rres = pow(werr*rwdot,2);
      double rpull = dr/rres;
      unsigned ist = hhit._flag.hasAnyProperty(StrawHitFlag::stereo) ? 0 : 1;
      // cut
      if( fabs(dr) > _maxdr || fabs(rpull) > _maxrpull || rwdot > _maxrwdot[ist]){
	// outlier hit flag it
	hhit._flag.merge(outlier);
      } else {
	// clear the hit in case it was formerly an outlier
	hhit._flag.clear(outlier);
      }
      changed |= oldout != hhit._flag.hasAnyProperty(outlier);
    }
    return changed;
  }

  bool RobustHelixFinder::filterHits(HelixSeed& hseed) {
  // 3d selection on top of radial selection
    RobustHelix& helix = hseed._helix;
    HelixHitCollection& hhits = hseed._hhits;
    bool changed(false);
    static StrawHitFlag outlier(StrawHitFlag::outlier);
    static Hep3Vector zaxis(0.0,0.0,1.0); // unit in z direction
// loop over hits
    for(auto& hhit : hhits) {
    // only remove hits, never add
      if(!hhit._flag.hasAnyProperty(outlier)) {
	// first compare detector azimuth with the circle center
	double hphi = hhit.pos().phi(); 
	double dphi = fabs(Angles::deltaPhi(hphi,helix.fcent()));
	// compute spatial distance and chisquiared
	// directions
	Hep3Vector const& wdir = hhit.wdir();
	Hep3Vector wtdir = zaxis.cross(wdir); // transverse direction to the wire
	Hep3Vector cvec = (hhit.pos() - helix.center()).perpPart(); // direction from the circle center to the hit
	Hep3Vector cdir = cvec.unit(); // direction from the circle center to the hit
	Hep3Vector cperp = zaxis.cross(cdir); // direction perp to the radius
	// positions
	Hep3Vector hpos = hhit.pos(); // this sets the z position to the hit z
	helix.position(hpos); // this computes the helix expectation at that z
	Hep3Vector dh = hhit.pos() - hpos; // this is the vector between them
	double dtrans = fabs(dh.dot(wtdir)); // transverse projection
	double dwire = fabs(dh.dot(wdir)); // projection along wire direction
	// compute the total resolution including hit and helix parameters first along the wire
	double wres2 = std::pow(hhit.posRes(StrawHitPosition::wire),(int)2) +
	  std::pow(_cradres*cdir.dot(wdir),(int)2) +
	  std::pow(_cperpres*cperp.dot(wdir),(int)2);
	// transverse to the wires
	double wtres2 = std::pow(hhit.posRes(StrawHitPosition::trans),(int)2) +
	  std::pow(_cradres*cdir.dot(wtdir),(int)2) +
	  std::pow(_cperpres*cperp.dot(wtdir),(int)2);
	// create a chisquared from these
	double chisq = sqrt( dwire*dwire/wres2 + dtrans*dtrans/wtres2 );
	// cut
	if( dphi > _maxphisep || fabs(dwire) > _maxdwire || fabs(dtrans) > _maxdtrans || chisq > _maxchisq) {
	  // outlier hit flag it.  But don't clear if it's not an outlier
	  hhit._flag.merge(outlier);
	  changed = true;
	}
      }
    }
    return changed;
  }

  bool RobustHelixFinder::findData(const art::Event& evt){
    _shcol = 0; _shfcol = 0; _shpcol = 0; _tccol = 0; 
    auto shH = evt.getValidHandle<StrawHitCollection>(_shTag);
    _shcol = shH.product();
    auto shpH = evt.getValidHandle<StrawHitPositionCollection>(_shpTag);
    _shpcol = shpH.product();
    auto shfH = evt.getValidHandle<StrawHitFlagCollection>(_shfTag);
    _shfcol = shfH.product();
    auto tcH = evt.getValidHandle<TimeClusterCollection>(_tcTag);
    _tccol = tcH.product();

    return _shcol != 0 && _shfcol != 0 && _shpcol != 0 && _tccol != 0;
  }

  void RobustHelixFinder::prefilterHits(HelixSeed& hseed ) {
    static StrawHitFlag outlier(StrawHitFlag::outlier);
    HelixHitCollection& hhits = hseed._hhits;
    // Iteratively use the average X and Y to define the average phi of the hits
    bool changed(true);
    size_t nhit = hhits.size();
    while(changed && nhit > 0){
      nhit = 0;
      changed = false;
      accumulator_set<double, stats<tag::median(with_p_square_quantile) > > accx;
      accumulator_set<double, stats<tag::median(with_p_square_quantile) > > accy;
      for(auto const& hhit : hhits ) {
	if(!hhit._flag.hasAnyProperty(outlier)){
	  accx(hhit._pos.x());
	  accy(hhit._pos.y());
	  nhit++;
	}
      }
      // compute median phi
      double mx = extract_result<tag::median>(accx);
      double my = extract_result<tag::median>(accy);
      double mphi = atan2(my,mx);
      // find the worst hit
      double maxdphi =0.0;
      auto worsthit = hhits.end();
      for(auto ihit = hhits.begin(); ihit != hhits.end(); ++ihit) {
	if(!ihit->_flag.hasAnyProperty(outlier)){
	  double phi = ihit->pos().phi();
	  double dphi = fabs(Angles::deltaPhi(phi,mphi));
	  if(dphi > maxdphi){
	    maxdphi = dphi;
	    worsthit = ihit;
	  }
	}
      }
      if(maxdphi > _maxphisep){
	worsthit->_flag.merge(outlier);
	changed = true;
      }
    }
  }

  void RobustHelixFinder::updateT0(HelixSeed& hseed) {
    static StrawHitFlag outlier(StrawHitFlag::outlier); // make this a class member FIXME!
  // this code should be stanardized and moved to TrkTimeCalculator FIXME!
    auto const& hhits = hseed.hits();
    accumulator_set<double, stats<tag::weighted_variance(lazy)>, double > terr;
    for(auto const& hhit : hhits) {
      if(!hhit._flag.hasAnyProperty(outlier)){
	// weight inversely by the hit time resolution
	double wt = std::pow(1.0/_ttcalc.strawHitTimeErr(),2);
	StrawHitIndex ind = hhit.index();
	terr(_ttcalc.strawHitTime(_shcol->at(ind),_shpcol->at(ind)),weight=wt);
      }
    }
    // add cluster time 
    if(hseed.caloCluster().isNonnull()){
      double time = _ttcalc.caloClusterTime(*hseed.caloCluster());
      double wt = std::pow(1.0/_ttcalc.caloClusterTimeErr(hseed.caloCluster()->sectionId()),2);
      terr(time,weight=wt);
    }
    hseed._t0._t0 = extract_result<tag::weighted_mean>(terr) + _t0shift; // ad-hoc correction FIXME!!
    hseed._t0._t0err = sqrt(std::max(0.0,extract_result<tag::weighted_variance(lazy)>(terr))/extract_result<tag::count>(terr));
  }

  void RobustHelixFinder::correctParameters(RobustHelix& helix) {
  // linear correction to center radius
    double cr = helix.rcent();
    double crcorr = _crcorr[0] + _crcorr[1]*cr;
  // quadratic correciton to radius (relative to a linear fit)
    double dr = helix.radius()-_rcorr[0];
    static double qf = pow(1.0+_rcorr[1]*_rcorr[1],1.5)*_rcorr[2];
    double rcorr = _rcorr[3] + _rcorr[1]*dr + qf*dr*dr;
    // update
    helix._rcent = crcorr;
    helix._radius = rcorr;
  }

  void RobustHelixFinder::fitHelix(HelixSeed& hseed) {
    // iteratively fit the circle
    unsigned niter(0);
    bool changed(true);
    do {
      _hfit.fitCircle(hseed);
      changed = filterCircleHits(hseed);
      ++niter;
    } while(hseed._status.hasAllProperties(TrkFitFlag::circleOK)  && niter < _maxniter && changed);
    // then fit phi-Z
    if(hseed._status.hasAnyProperty(TrkFitFlag::circleOK)){
      if(niter < _maxniter) hseed._status.merge(TrkFitFlag::circleConverged);
      // solve for the longitudinal parameters
      niter = 0;
      changed = false;
      do {
	_hfit.fitFZ(hseed);
	changed = filterHits(hseed);
	++niter;
      } while(hseed._status.hasAllProperties(TrkFitFlag::phizOK)  && niter < _maxniter && changed);
      if(hseed._status.hasAnyProperty(TrkFitFlag::phizOK)){
	if(niter < _maxniter) hseed._status.merge(TrkFitFlag::phizConverged);
	// final test
	if (_hfit.goodHelix(hseed.helix()))hseed._status.merge(TrkFitFlag::helixOK);
	if(_diag > 0)_niter->Fill(niter);
      }
    }
  }

  unsigned RobustHelixFinder::hitCount(HelixSeed const& hseed) {
    static StrawHitFlag outlier(StrawHitFlag::outlier);
    HelixHitCollection const& hhits = hseed._hhits;
    unsigned retval(0);
    for(auto hhit : hhits)
      if(!hhit._flag.hasAnyProperty(outlier))++retval;
    return retval;
  }

}
using mu2e::RobustHelixFinder;
DEFINE_ART_MODULE(RobustHelixFinder);
