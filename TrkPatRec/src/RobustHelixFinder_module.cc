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
    Double_t& _drho; // hit transverse radius minus helix radius
    Double_t& _dphi; // hit azimuth minus helix azimuth (at the hit z)
    Double_t& _hrho; // helix transverse radius (at the hit z)
    Double_t& _chisq; // chisq of spatial information, using average errors
    Double_t& _dt;  // time difference of hit WRT average
    Double_t& _whdot; // cosine of the angle between the wire and the helix tangent
    Double_t& _hhrho; // hit transverse radius
    HelixHitMVA() : _pars(9,0.0),_dtrans(_pars[0]),_dwire(_pars[1]),_drho(_pars[2]),_dphi(_pars[3]),_hrho(_pars[4]),_chisq(_pars[5]),_dt(_pars[6]),_whdot(_pars[7]),_hhrho(_pars[8]){}
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
    double				_maxphisep; // maximum separation in global azimuth of hits
    bool				_saveall; // write out all helices regardless of fit success
    unsigned				_maxniter;  // maximum # of iterations over outlier filtering + fitting
    double				_cradres; // average center resolution along center position (mm)
    double				_cperpres; // average center resolution perp to center position (mm)
    double				_radres; // average radial resolution for circle fit (mm)
    double				_phires; // average azimuthal resolution on circle (rad)
    double				_maxdwire; // outlier cut on distance between hit and helix along wire
    double				_maxdtrans; // outlier cut on distance between hit and helix perp to wire
    double				_maxchisq; // outlier cut on chisquared
    bool				_usemva; // use MVA or not.  If not, simple cuts
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

    // cache of event objects
    const StrawHitCollection*  _shcol;
    const StrawHitPositionCollection*  _shpcol;
    const StrawHitFlagCollection*      _shfcol;
    const TimeClusterCollection*       _tccol;

    // robust helix fitter
    RobustHelixFit                     _hfit;
    TrkTimeCalculator			_ttcalc;
    bool			      _parcor; // correct the parameters
    double _rccorr[2], _rcorr[2]; // corrections for center radius and radius

    // helper functions
    bool findData           (const art::Event& e);
    void prefilterHits(HelixSeed& hseed); // rough filtering based on hits being in 1/2 the detector
    bool filterHits(HelixSeed& hseed); // return value tells if any hits changed state
    void updateT0(HelixSeed& hseed); // update T0 value based on current good hits
    void correctParameters(RobustHelix& helix);// correct the helix parameters to MC truth using linear relations and correlations
 
  };

  RobustHelixFinder::RobustHelixFinder(fhicl::ParameterSet const& pset) :
    _diag        (pset.get<int>("diagLevel",0)),
    _debug       (pset.get<int>("debugLevel",0)),
    _printfreq   (pset.get<int>("printFrequency",101)),
    _prefilter   (pset.get<bool>("PrefilterHits",true)),
    _maxphisep(pset.get<double>("MaxPhiHitSeparation",1.0)),
    _saveall     (pset.get<bool>("SaveAllHelices",false)),
    _maxniter    (pset.get<unsigned>("MaxIterations",0)), // iterations over outlier removal
    _cradres	 (pset.get<double>("CenterRadialResolution",12.0)),
    _cperpres	 (pset.get<double>("CenterPerpResolution",12.0)),
    _radres	 (pset.get<double>("RadiusResolution",10.0)),
    _phires	 (pset.get<double>("AzimuthREsolution",0.1)),
    _maxdwire    (pset.get<double>("MaxWireDistance",200.0)), // max distance along wire
    _maxdtrans   (pset.get<double>("MaxTransDistance",100.0)), // max distance perp to wire (and z)
    _maxchisq    (pset.get<double>("MaxChisquared",9.0)), // max chisquared
    _usemva      (pset.get<bool>("UseMVAHitSelection",true)),
    _minmva      (pset.get<double>("MinMVA",0.05)), // min MVA output
    _shTag	 (pset.get<art::InputTag>("StrawHitCollection","makeSH")),
    _shpTag	 (pset.get<art::InputTag>("StrawHitPositionCollection","MakeStereoHits")),
    _shfTag	 (pset.get<art::InputTag>("StrawHitFlagCollection","TimeClusterFinder")),
    _tcTag	 (pset.get<art::InputTag>("TimeClusterCollection","TimeClusterFinder")),
    _trackseed   (pset.get<string>("HelixSeedCollectionLabel","TimeClusterFinder")),
    _hsel        (pset.get<std::vector<std::string> >("HitSelectionBits")),
    _hbkg        (pset.get<vector<string> >("HitBackgroundBits",vector<string>{"DeltaRay","Isolated"})),
    _mvatool     (pset.get<fhicl::ParameterSet>("HelixHitMVA",fhicl::ParameterSet())),
    _hfit        (pset.get<fhicl::ParameterSet>("RobustHelixFit",fhicl::ParameterSet())),
    _ttcalc            (pset.get<fhicl::ParameterSet>("T0Calculator",fhicl::ParameterSet())),
    _parcor      (pset.get<bool>("CorrectHelixParameters",false))
  {
    _rccorr[0] = pset.get<double>("CenterRadiusOffset",118.0); // mm
    _rccorr[1] = pset.get<double>("CenterRadiusSlope",0.57); // RC_reco = _rccorr[0] + _rccorr[1]*RC_MC
    _rcorr[0] = pset.get<double>("RadiusOffset",74.0); // mm
    _rcorr[1] = pset.get<double>("RadiusSlope",0.73); // R_reco = _rcorr[0] + _rcorr[1]*R_MC

    produces<HelixSeedCollection>();
  }

  RobustHelixFinder::~RobustHelixFinder(){}

  void RobustHelixFinder::beginJob() {
  // initialize MVA
    if(_usemva) {
      _mvatool.initMVA();
      if(_debug > 0){
	cout << "RobustHeilxFinder MVA parameters: " << endl;
	_mvatool.showMVA();
      }
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
      // iteratively fit the helix and filter outliers in space and time
      unsigned niter(0);
      bool changed(true);
      do {
	++niter;
	_hfit.fitHelix(hseed);
	changed = filterHits(hseed);
	if(changed)updateT0(hseed);
      } while(hseed._status.hasAllProperties(TrkFitFlag::helixOK)  && niter < _maxniter && changed);
      if(niter < _maxniter)hseed._status.merge(TrkFitFlag::helixConverged);
      // final test
      if((hseed._status.hasAllProperties(TrkFitFlag::helixOK) && _hfit.helicity() == hseed._helix.helicity()) || _saveall) {
// correct the parameters if requested
	if(_parcor)correctParameters(hseed._helix);
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

  bool RobustHelixFinder::filterHits(HelixSeed& hseed) {
    RobustHelix& helix = hseed._helix;
    HelixHitCollection& hhits = hseed._hhits;
    bool changed(false);
    static StrawHitFlag outlier(StrawHitFlag::outlier);
    static Hep3Vector zaxis(0.0,0.0,1.0); // unit in z direction
// loop over hits
    for(auto& hhit : hhits) {
      bool oldout = !hhit._flag.hasAnyProperty(outlier);
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
      // fill MVA struct
      _vmva._dtrans = dh.dot(wtdir); // transverse projection
      _vmva._dwire = dh.dot(wdir); // projection along wire direction
      _vmva._drho = cvec.mag() - helix.radius(); // radius difference
      _vmva._whdot = wdir.dot(cperp); // cosine of angle between wire and helix tangent 
      _vmva._dphi = hhit._phi - helix.circleAzimuth(hhit.pos().z()); // azimuth difference WRT circle center
      _vmva._hhrho = cvec.mag(); // hit transverse radius WRT circle center
      _vmva._hrho = hpos.perp();  // hit detector transverse radius
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
      // compute the MVA
      double mvaout = _mvatool.evalMVA(_vmva._pars);
      // update the hit
      hhit._hqual = mvaout;
      if( dphi > _maxphisep || fabs(_vmva._dwire) > _maxdwire || fabs(_vmva._dtrans) > _maxdtrans || _vmva._chisq > _maxchisq || mvaout < _minmva) {
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
    hseed._t0._t0 = extract_result<tag::weighted_mean>(terr);
    hseed._t0._t0err = sqrt(std::max(0.0,extract_result<tag::weighted_variance(lazy)>(terr))/extract_result<tag::count>(terr));
  }

  void RobustHelixFinder::correctParameters(RobustHelix& helix) {
    // invert the parameters, since we want the best match with truth
    helix._radius = (helix._radius - _rcorr[0])/_rcorr[1];
    helix._rcent = (helix._rcent - _rccorr[0])/_rccorr[1];
  }
}
using mu2e::RobustHelixFinder;
DEFINE_ART_MODULE(RobustHelixFinder);
