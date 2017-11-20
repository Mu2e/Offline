//
// TTracker Pattern Recognition based on Robust Helix Fit
//
// $Id: RobustHelixFinder_module.cc,v 1.2 2014/08/30 12:19:38 tassiell Exp $
// $Author: tassiell $
// $Date: 2014/08/30 12:19:38 $
//
// Original author D. Brown and G. Tassielli
//

#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "GeneralUtilities/inc/Angles.hh"
#include "Mu2eUtilities/inc/MVATools.hh"

#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/StereoHit.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
#include "RecoDataProducts/inc/HelixSeed.hh"
#include "RecoDataProducts/inc/TrkFitFlag.hh"

#include "TrkReco/inc/TrkTimeCalculator.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "BTrk/BaBar/BaBar.hh"
#include "TrkReco/inc/TrkDef.hh"
#include "BTrkData/inc/TrkStrawHit.hh"
#include "TrkReco/inc/RobustHelixFit.hh"

#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"

#include <boost/accumulators/accumulators.hpp>
#include "boost_fix/accumulators/statistics/stats.hpp"
#include "boost_fix/accumulators/statistics.hpp"
#include <boost/accumulators/statistics/median.hpp>

#include "TH1F.h"

#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <utility>
#include <functional>
#include <float.h>
#include <vector>
#include <map>

using namespace std;
using namespace boost::accumulators;


namespace {

  struct HelixHitMVA
  {
    std::vector <Double_t> _pars;
    Double_t& _dtrans; // distance from hit to helix perp to the wrire
    Double_t& _dwire;  // distance from hit to helix along the wrire
    Double_t& _chisq;  // chisq of spatial information, using average errors
    Double_t& _dt;     // time difference of hit WRT average
    Double_t& _drho;   // hit transverse radius minus helix radius
    Double_t& _dphi;   // hit azimuth minus helix azimuth (at the hit z)
    Double_t& _rwdot;  // dot product between circle radial direction and wire direction
    Double_t& _hrho;   // helix transverse radius (at the hit z)
    Double_t& _hhrho;  // hit transverse radius
    HelixHitMVA() : _pars(9,0.0),_dtrans(_pars[0]),_dwire(_pars[1]),_chisq(_pars[2]),_dt(_pars[3]),
     _drho(_pars[4]),_dphi(_pars[5]),_rwdot(_pars[6]),_hrho(_pars[7]),_hhrho(_pars[8]){}
  };

}

namespace mu2e {

  class RobustHelixFinder : public art::EDProducer {
  public:
    explicit RobustHelixFinder(fhicl::ParameterSet const&);
    virtual ~RobustHelixFinder();
    virtual void beginJob();
    virtual void produce(art::Event& event );

  private:
    int                                 _diag,_debug;
    int                                 _printfreq;
    bool				_prefilter; // prefilter hits based on sector
    bool				_updatestereo; // update the stereo hit positions each iteration
    double				_dhit; // distance between hit position updates to consider changed' (mm)
    double				_dhit2;
    unsigned				_minnhit; // minimum # of hits to work with
    double				_maxdr; // maximum hit-helix radius difference
    double				_maxrpull; // maximum hit-helix radius difference pull
    double				_maxphisep; // maximum separation in global azimuth of hits
    TrkFitFlag				_saveflag; // write out all helices that satisfy these flags
    unsigned				_maxniter;  // maximum # of iterations over outlier filtering + fitting
    double				_cradres; // average center resolution along center position (mm)
    double				_cperpres; // average center resolution perp to center position (mm)
    double				_maxdwire; // outlier cut on distance between hit and helix along wire
    double				_maxdtrans; // outlier cut on distance between hit and helix perp to wire
    double				_maxchisq; // outlier cut on chisquared
    double				_maxrwdot[2]; // outlier cut on angle between radial direction and wire: smaller is better
    double				_minrerr; // minimum radius error

    bool				_usemva; // use MVA to cut outliers
    double				_minmva; // outlier cut on MVA

    art::InputTag			_ccFastTag;
    art::InputTag			_shTag;
    art::InputTag			_shpTag;
    art::InputTag			_shfTag;
    art::InputTag			_tcTag;
    art::InputTag			_sthTag;
    std::string                         _trackseed;

    StrawHitFlag  _hsel, _hbkg;

    MVATools _stmva, _nsmva;
    HelixHitMVA _vmva; // input variables to TMVA for filtering hits

    TH1F* _niter, *_nitermva;

    RobustHelixFit   _hfit;
    std::vector<Helicity> _hels; // helicity values to fit 
    TrkTimeCalculator _ttcalc;
    double            _t0shift;   
    StrawHitFlag      _outlier;
    bool              _updateStereo;
    
    
    void     findHelices(const StrawHitCollection& shcol,const StrawHitPositionCollection& shpcol, 
                         const StrawHitFlagCollection& shfcol,const StereoHitCollection& sthcol, const TimeClusterCollection& tccol, 
                         HelixSeedCollection& outseeds, art::Event& event);    
    void     prefilterHits(HelixSeed& hseed); 
    unsigned filterCircleHits(HelixSeed& hseed); 
    bool     filterHits(HelixSeed& hseed);
    void     fillMVA(HelixSeed& hseed, const StrawHitCollection& shcol); 
    bool     filterHitsMVA(HelixSeed& hseed);
    void     updateT0(HelixSeed& hseed, const StrawHitCollection& shcol,const StrawHitPositionCollection& shpcol);
    bool     updateStereo(HelixSeed& hseed, const StereoHitCollection& sthcol, const StrawHitCollection& shcol);
    unsigned hitCount(HelixSeed const& hseed);
    void     fitHelix(HelixSeed& hseed, const StrawHitCollection& shcol,const StrawHitPositionCollection& shpcol, const StereoHitCollection& sthcol);
    void     refitHelix(HelixSeed& hseed);

  };

  RobustHelixFinder::RobustHelixFinder(fhicl::ParameterSet const& pset) :
    _diag        (pset.get<int>("diagLevel",0)),
    _debug       (pset.get<int>("debugLevel",0)),
    _printfreq   (pset.get<int>("printFrequency",101)),
    _prefilter   (pset.get<bool>("PrefilterHits",true)),
    _updatestereo(pset.get<bool>("UpdateStereoHits",true)),
    _dhit	 (pset.get<double>("HitDistanceChange",10.0)), // mm
    _minnhit	 (pset.get<unsigned>("minNHit",5)),
    _maxdr	 (pset.get<double>("MaxRadiusDiff",100.0)), // mm
    _maxrpull	 (pset.get<double>("MaxRPull",5.0)), // unitless
    _maxphisep	 (pset.get<double>("MaxPhiHitSeparation",1.0)),
    _saveflag    (pset.get<vector<string> >("SaveHelixFlag",vector<string>{"HelixOK"})),
    _maxniter    (pset.get<unsigned>("MaxIterations",10)), // iterations over outlier removal
    _cradres	 (pset.get<double>("CenterRadialResolution",20.0)),
    _cperpres	 (pset.get<double>("CenterPerpResolution",12.0)),
    _maxdwire    (pset.get<double>("MaxWireDistance",200.0)), // max distance along wire
    _maxdtrans   (pset.get<double>("MaxTransDistance",80.0)), // max distance perp to wire (and z)
    _maxchisq    (pset.get<double>("MaxChisquared",100.0)), // max chisquared
    _minrerr     (pset.get<double>("MinRadiusErr",20.0)), // mm
    _usemva      (pset.get<bool>("UseHitMVA",true)),
    _minmva      (pset.get<double> ("MinMVA",0.1)), // min MVA output to define an outlier
    _ccFastTag	 (pset.get<art::InputTag>("CaloCluster","CaloClusterFast")),
    _shTag	 (pset.get<art::InputTag>("StrawHitCollection","makeSH")),
    _shpTag	 (pset.get<art::InputTag>("StrawHitPositionCollection","MakeStereoHits")),
    _shfTag	 (pset.get<art::InputTag>("StrawHitFlagCollection","TimeClusterFinder")),
    _tcTag	 (pset.get<art::InputTag>("TimeClusterCollection","TimeClusterFinder")),
    _sthTag	 (pset.get<art::InputTag>("StereoHitCollection","MakeStereoHits")),
    _trackseed   (pset.get<string>("HelixSeedCollectionLabel","TimeClusterFinder")),
    _hsel        (pset.get<std::vector<std::string> >("HitSelectionBits")),
    _hbkg        (pset.get<std::vector<std::string> >("HitBackgroundBits",std::vector<std::string>{"Background"})),
    _stmva       (pset.get<fhicl::ParameterSet>("HelixStereoHitMVA",fhicl::ParameterSet())),
    _nsmva       (pset.get<fhicl::ParameterSet>("HelixNonStereoHitMVA",fhicl::ParameterSet())),
    _hfit        (pset.get<fhicl::ParameterSet>("RobustHelixFit",fhicl::ParameterSet())),
    _ttcalc      (pset.get<fhicl::ParameterSet>("T0Calculator",fhicl::ParameterSet())),
    _t0shift     (pset.get<double>("T0Shift",4.0)),
    _outlier     (StrawHitFlag::outlier),
    _updateStereo    (pset.get<bool>("UpdateStereo",true))
  {
    _maxrwdot[0] = pset.get<double>("MaxStereoRWDot",1.0);
    _maxrwdot[1] = pset.get<double>("MaxNonStereoRWDot",1.0);
    _dhit2 = _dhit*_dhit;
    std::vector<int> helvals = pset.get<std::vector<int> >("Helicities",vector<int>{Helicity::neghel,Helicity::poshel});
    for(auto hv : helvals) {
      Helicity hel(hv);
      _hels.push_back(hel);
      produces<HelixSeedCollection>(Helicity::name(hel));
    }
  }

  RobustHelixFinder::~RobustHelixFinder(){}

  void RobustHelixFinder::beginJob() {
    _stmva.initMVA();
    _nsmva.initMVA();
    if (_debug > 0)
    {
      std::cout << "RobustHeilxFinder Stereo Hit MVA parameters: " << std::endl;
      _stmva.showMVA();
      std::cout << "RobustHeilxFinder Non-Stereo Hit MVA parameters: " << std::endl;
      _nsmva.showMVA();
    }

    if (_diag > 0){
      art::ServiceHandle<art::TFileService> tfs;
      _niter = tfs->make<TH1F>( "niter" , "Number of Fit Iteraions",201,-0.5,200.5);
      _nitermva = tfs->make<TH1F>( "nitermva" , "Number of MVA Fit Iteraions",201,-0.5,200.5);
    }
  }

  void RobustHelixFinder::produce(art::Event& event ) {
    // find input
    auto tcH = event.getValidHandle<TimeClusterCollection>(_tcTag);
    const TimeClusterCollection& tccol(*tcH);

    auto shH = event.getValidHandle<StrawHitCollection>(_shTag);
    const StrawHitCollection& shcol(*shH);

    auto spH = event.getValidHandle<StrawHitPositionCollection>(_shpTag);
    const StrawHitPositionCollection& shpcol(*spH);

    auto sfH = event.getValidHandle<StrawHitFlagCollection>(_shfTag);
    const StrawHitFlagCollection& shfcol(*sfH);

    auto stH = event.getValidHandle<StereoHitCollection>(_sthTag);
    const StereoHitCollection& sthcol(*stH);

    // create output: seperate by helicity
    std::map<Helicity,unique_ptr<HelixSeedCollection>> helcols;
    for( auto const& hel : _hels) 
      helcols[hel] = unique_ptr<HelixSeedCollection>(new HelixSeedCollection());

    // create initial helicies from time clusters: to begin, don't specificy helicity
    for (size_t index=0;index< tccol.size();++index) {
      const auto& tclust = tccol.at(index);
      HelixSeed hseed;
      hseed._t0 = tclust._t0;
      hseed._timeCluster = art::Ptr<TimeCluster>(tcH,index);
// convert straw hits to helix hits: these cache additional information needed to improve efficiency
      for (const auto& ind : tclust._strawHitIdxs) {
	if (shfcol.at(ind).hasAnyProperty(_hsel) && !shfcol.at(ind).hasAnyProperty(_hbkg)) {
	  HelixHit hhit(shpcol.at(ind),ind);
	  hhit._flag.clear(StrawHitFlag::resolvedphi);
	  hhit._flag.merge(shfcol.at(ind));
	  hseed._hhits.push_back(std::move(hhit));        
	}  
      }
      // filter hits and test
      if (_prefilter) prefilterHits(hseed);
      if (hitCount(hseed) >= _minnhit){
	hseed._status.merge(TrkFitFlag::hitsOK);
	// initial circle fit
	_hfit.fitCircle(hseed);
	if (hseed._status.hasAnyProperty(TrkFitFlag::circleOK)) {
	  // loop over helicities. 
	  for(auto const& hel : _hels ) {
	  // tentatively put a copy with the specified helicity in the appropriate output vector
	    hseed._helix._helicity = hel;
	    HelixSeedCollection* hcol = helcols[hel].get();
	    hcol->push_back(hseed);
	  // attempt complete fit 
	    fitHelix(hcol->back(),shcol, shpcol, sthcol);
	    // test fit status; if not successful, pop it off
	    if (!hcol->back().status().hasAllProperties(_saveflag))
	      hcol->pop_back();
	  }
	}	
      }	
    }
    // put final collections into event 
    for(auto const& hel : _hels ) {
      event.put(std::move(helcols[hel]),Helicity::name(hel));
    }
  }

  void RobustHelixFinder::fillMVA(HelixSeed& hseed, const StrawHitCollection& shcol)
  {
    RobustHelix& helix = hseed._helix;

    static CLHEP::Hep3Vector zaxis(0.0,0.0,1.0); // unit in z direction
    for(auto& hhit : hseed._hhits)
    {

      const CLHEP::Hep3Vector& wdir = hhit.wdir();
      CLHEP::Hep3Vector wtdir = zaxis.cross(wdir); // transverse direction to the wire
      CLHEP::Hep3Vector cvec = (hhit.pos() - helix.center()).perpPart(); // direction from the circle center to the hit
      CLHEP::Hep3Vector cdir = cvec.unit();        // direction from the circle center to the hit
      CLHEP::Hep3Vector cperp = zaxis.cross(cdir); // direction perp to the radius

      CLHEP::Hep3Vector hpos = hhit.pos();      // this sets the z position to the hit z
      helix.position(hpos);                     // this computes the helix expectation at that z
      CLHEP::Hep3Vector dh = hhit.pos() - hpos; // this is the vector between them

      _vmva._dtrans = fabs(dh.dot(wtdir));              // transverse projection
      _vmva._dwire = fabs(dh.dot(wdir));               // projection along wire direction
      _vmva._drho = fabs(cvec.mag() - helix.radius()); // radius difference
      _vmva._dphi = fabs(hhit._phi - helix.circleAzimuth(hhit.pos().z())); // azimuth difference WRT circle center
      _vmva._hhrho = cvec.mag();            // hit transverse radius WRT circle center
      _vmva._hrho = hpos.perp();            // hit detector transverse radius
      _vmva._rwdot = fabs(wdir.dot(cdir));  // compare directions of radius and wire

      // compute the total resolution including hit and helix parameters first along the wire
      double wres2 = std::pow(hhit.posRes(StrawHitPosition::wire),(int)2) +
	std::pow(_cradres*cdir.dot(wdir),(int)2) +
	std::pow(_cperpres*cperp.dot(wdir),(int)2);

      // transverse to the wires
      double wtres2 = std::pow(hhit.posRes(StrawHitPosition::trans),(int)2) +
	std::pow(_cradres*cdir.dot(wtdir),(int)2) +
	std::pow(_cperpres*cperp.dot(wtdir),(int)2);

      _vmva._chisq = sqrt( _vmva._dwire*_vmva._dwire/wres2 + _vmva._dtrans*_vmva._dtrans/wtres2 );          
      _vmva._dt = shcol.at(hhit._shidx).time() - hseed._t0.t0();

      if (hhit._flag.hasAnyProperty(StrawHitFlag::stereo))
      {
	hhit._hqual = _stmva.evalMVA(_vmva._pars);
      } else {
	hhit._hqual = _nsmva.evalMVA(_vmva._pars);
      }
    }
  }

  bool RobustHelixFinder::filterHitsMVA(HelixSeed& hseed)
  {  
    bool changed(false);
    for (auto& hhit : hseed._hhits)
    {
      bool oldout = hhit._flag.hasAnyProperty(_outlier);

      if (hhit._hqual < _minmva ) hhit._flag.merge(_outlier);
      else                        hhit._flag.clear(_outlier);

      changed |= oldout != hhit._flag.hasAnyProperty(_outlier);
    }
    return changed;
  }

  unsigned  RobustHelixFinder::filterCircleHits(HelixSeed& hseed)
  {
    unsigned changed(0);
    static CLHEP::Hep3Vector zaxis(0.0,0.0,1.0); // unit in z direction
    RobustHelix& helix = hseed._helix;

    for(auto& hhit : hseed._hhits)
    {
      bool oldout = hhit._flag.hasAnyProperty(_outlier);

      const CLHEP::Hep3Vector& wdir = hhit.wdir();
      CLHEP::Hep3Vector cvec = (hhit.pos() - helix.center()).perpPart(); // direction from the circle center to the hit
      CLHEP::Hep3Vector cdir = cvec.unit(); // direction from the circle center to the hit
      double rwdot = wdir.dot(cdir); // compare directions of radius and wire
      double rwdot2 = rwdot*rwdot;
      // compute radial difference and pull
      double dr = cvec.mag()-helix.radius();
      double werr = hhit.posRes(StrawHitPosition::wire);
      double terr = hhit.posRes(StrawHitPosition::trans);
      // the resolution is dominated the resolution along the wire
      double rres = std::max(sqrt(werr*werr*rwdot2 + terr*terr*(1.0-rwdot2)),_minrerr);
      double rpull = dr/rres;
      unsigned ist = hhit._flag.hasAnyProperty(StrawHitFlag::stereo) ? 0 : 1;

      if ( std::abs(dr) > _maxdr || fabs(rpull) > _maxrpull || rwdot > _maxrwdot[ist])
      {
	hhit._flag.merge(_outlier);
      } else {
	hhit._flag.clear(_outlier);
      }
      if (oldout != hhit._flag.hasAnyProperty(_outlier)) ++changed;
    }
    return changed;
  }


  // 3d selection on top of radial selection
  bool RobustHelixFinder::filterHits(HelixSeed& hseed)
  {
    RobustHelix& helix = hseed._helix;
    HelixHitCollection& hhits = hseed._hhits;
    bool changed(false);
    static CLHEP::Hep3Vector zaxis(0.0,0.0,1.0); // unit in z direction

    // loop over hits
    for(auto& hhit : hhits)
    {     
      if (hhit._flag.hasAnyProperty(_outlier)) continue;

      double hphi = hhit.pos().phi();
      double dphi = fabs(Angles::deltaPhi(hphi,helix.fcent()));

      const CLHEP::Hep3Vector& wdir = hhit.wdir();
      CLHEP::Hep3Vector wtdir = zaxis.cross(wdir);   // transverse direction to the wire
      CLHEP::Hep3Vector cvec = (hhit.pos() - helix.center()).perpPart(); // direction from the circle center to the hit
      CLHEP::Hep3Vector cdir = cvec.unit();          // direction from the circle center to the hit
      CLHEP::Hep3Vector cperp = zaxis.cross(cdir);   // direction perp to the radius

      CLHEP::Hep3Vector hpos = hhit.pos(); // this sets the z position to the hit z
      helix.position(hpos);                // this computes the helix expectation at that z
      CLHEP::Hep3Vector dh = hhit.pos() - hpos;   // this is the vector between them
      double dtrans = fabs(dh.dot(wtdir)); // transverse projection
      double dwire = fabs(dh.dot(wdir));   // projection along wire direction

      // compute the total resolution including hit and helix parameters first along the wire
      double wres2 = std::pow(hhit.posRes(StrawHitPosition::wire),(int)2) +
	std::pow(_cradres*cdir.dot(wdir),(int)2) +
	std::pow(_cperpres*cperp.dot(wdir),(int)2);
      // transverse to the wires
      double wtres2 = std::pow(hhit.posRes(StrawHitPosition::trans),(int)2) +
	std::pow(_cradres*cdir.dot(wtdir),(int)2) +
	std::pow(_cperpres*cperp.dot(wtdir),(int)2);

      double chisq = dwire*dwire/wres2 + dtrans*dtrans/wtres2;

      if( dphi > _maxphisep || fabs(dwire) > _maxdwire || fabs(dtrans) > _maxdtrans || chisq > _maxchisq) 
      {
	hhit._flag.merge(_outlier);
	changed = true;
      }
    }
    return changed;
  }

  void RobustHelixFinder::prefilterHits(HelixSeed& hseed)
  {
    HelixHitCollection& hhits = hseed._hhits;

    bool changed(true);
    size_t nhit = hhits.size();
    while (changed && nhit > 0)
    {
      nhit = 0;
      changed = false;
      accumulator_set<double, stats<tag::median(with_p_square_quantile) > > accx;
      accumulator_set<double, stats<tag::median(with_p_square_quantile) > > accy;
      for (const auto& hhit : hhits ) 
      {
	if (hhit._flag.hasAnyProperty(_outlier)) continue;
	accx(hhit._pos.x());
	accy(hhit._pos.y());
	++nhit;
      }

      double mx = extract_result<tag::median>(accx);
      double my = extract_result<tag::median>(accy);
      double mphi = atan2(my,mx);

      double maxdphi{0.0};
      auto worsthit = hhits.end();
      for(auto ihit = hhits.begin(); ihit != hhits.end(); ++ihit)
      {
	if (ihit->_flag.hasAnyProperty(_outlier)) continue;
	double phi  = ihit->pos().phi();
	double dphi = fabs(Angles::deltaPhi(phi,mphi));
	if(dphi > maxdphi)
	{
	  maxdphi = dphi;
	  worsthit = ihit;
	}
      }

      if (maxdphi > _maxphisep)
      {
	worsthit->_flag.merge(_outlier);
	changed = true;
      }
    }
  }

  void RobustHelixFinder::updateT0(HelixSeed& hseed, const StrawHitCollection& shcol, const StrawHitPositionCollection& shpcol)
  {
    const auto& hhits = hseed.hits();
    accumulator_set<double, stats<tag::weighted_variance(lazy)>, double > terr;
    for (const auto& hhit : hhits) 
    {
      if (hhit._flag.hasAnyProperty(_outlier)) continue;
      double wt = std::pow(1.0/_ttcalc.strawHitTimeErr(),2);
      StrawHitIndex ind = hhit.index();
      terr(_ttcalc.strawHitTime(shcol.at(ind),shpcol.at(ind)),weight=wt);
    }

    if (hseed.caloCluster().isNonnull())
    {
      double time = _ttcalc.caloClusterTime(*hseed.caloCluster());
      double wt = std::pow(1.0/_ttcalc.caloClusterTimeErr(hseed.caloCluster()->diskId()),2);
      terr(time,weight=wt);
    }

    if (sum_of_weights(terr) > 0.0)
    {
      hseed._t0._t0 = extract_result<tag::weighted_mean>(terr) + _t0shift; // ad-hoc correction FIXME!!
      hseed._t0._t0err = sqrt(std::max(0.0,extract_result<tag::weighted_variance(lazy)>(terr))/extract_result<tag::count>(terr));
    }
  }

  void RobustHelixFinder::fitHelix(HelixSeed& hseed, const StrawHitCollection& shcol,const StrawHitPositionCollection& shpcol, const StereoHitCollection& sthcol){
    // iteratively fit the helix including filtering
    unsigned niter(0);
    unsigned nitermva(0);
    bool changed(true);

    do {
      unsigned xyniter(0);
      unsigned xychanged = filterCircleHits(hseed);
      while (hseed._status.hasAllProperties(TrkFitFlag::circleOK) && xyniter < _maxniter && xychanged) {
	_hfit.fitCircle(hseed);
	xychanged = filterCircleHits(hseed);
	++xyniter;
      } 
      // then fit phi-Z
      if (hseed._status.hasAnyProperty(TrkFitFlag::circleOK)) {
	if (xyniter < _maxniter)
	  hseed._status.merge(TrkFitFlag::circleConverged);
	else
	  hseed._status.clear(TrkFitFlag::circleConverged);

	// solve for the longitudinal parameters
	unsigned fzniter(0);
	bool fzchanged(false);
	do {
	  _hfit.fitFZ(hseed);
	  fzchanged = filterHits(hseed);
	  ++fzniter;
	} while (hseed._status.hasAllProperties(TrkFitFlag::phizOK)  && fzniter < _maxniter && fzchanged);

	if (hseed._status.hasAnyProperty(TrkFitFlag::phizOK)) {
	  if (fzniter < _maxniter)
	    hseed._status.merge(TrkFitFlag::phizConverged);
	  else
	    hseed._status.clear(TrkFitFlag::phizConverged);
	}
      }
      ++niter;

      // update the stereo hit positions; this checks how much the positions changed
      // do this only in non trigger mode

      if (_updateStereo && _hfit.goodHelix(hseed.helix())) changed = updateStereo(hseed, sthcol, shcol);      
    } while (_hfit.goodHelix(hseed.helix()) && niter < _maxniter && changed);


    if (_hfit.goodHelix(hseed.helix())) {
      hseed._status.merge(TrkFitFlag::helixOK);
      updateT0(hseed,shcol,shpcol);
      if (niter < _maxniter) hseed._status.merge(TrkFitFlag::helixConverged);

      if (_usemva) {
	bool changed = true;
	while (hseed._status.hasAllProperties(TrkFitFlag::helixOK)  && nitermva < _maxniter && changed) {
	  fillMVA(hseed,shcol);
	  changed = filterHitsMVA(hseed);
	  if (!changed) break;
	  refitHelix(hseed);
	  // update t0 each iteration as that's used in the MVA
	  updateT0(hseed,shcol,shpcol);
	  ++nitermva;
	}
	if (nitermva < _maxniter)
	  hseed._status.merge(TrkFitFlag::helixConverged);
	else
	  hseed._status.clear(TrkFitFlag::helixConverged);
      }
    }
    if (_diag > 0){
      _niter->Fill(niter);
      _nitermva->Fill(nitermva);
      if (!_usemva) fillMVA(hseed,shcol);
    }
  }

  void RobustHelixFinder::refitHelix(HelixSeed& hseed) {
    // reset the fit status flags, in case this is called iteratively
    hseed._status.clear(TrkFitFlag::helixOK);      
    _hfit.fitCircle(hseed);
    if (hseed._status.hasAnyProperty(TrkFitFlag::circleOK)) {
      _hfit.fitFZ(hseed);
      if (_hfit.goodHelix(hseed._helix)) hseed._status.merge(TrkFitFlag::helixOK);
    }
  }

  unsigned RobustHelixFinder::hitCount(const HelixSeed& hseed) {
    return std::count_if(hseed._hhits.begin(),hseed._hhits.end(),
	[&](const HelixHit& hhit){return !hhit.flag().hasAnyProperty(_outlier);});
  }


  bool RobustHelixFinder::updateStereo(HelixSeed& hseed, const StereoHitCollection& sthcol, const StrawHitCollection& shcol) {
    static StrawHitFlag stereo(StrawHitFlag::stereo);
    bool retval(false);

    if (_updatestereo) {
      const Tracker& tracker = getTrackerOrThrow();
      // create a map from stereo hits back to helix hits
      std::map<size_t, std::pair<int,int> > sthmap;
      for(size_t ihh = 0; ihh < hseed.hits().size(); ++ihh) {
	const auto& hhit = hseed.hits().at(ihh);
	if (hhit.flag().hasAnyProperty(stereo) && hhit.stereoHitIndex() >= 0) {
	  size_t stindex = hhit.stereoHitIndex();
	  const auto& sthit = sthcol.at(stindex);
	  // find out which hit this is WRT the stereo indexing (first or second)
	  bool isfirst;
	  if (sthit.hitIndex1() == hhit.index())
	    isfirst = true;
	  else if (sthit.hitIndex2() == hhit.index())
	    isfirst = false;
	  else
	    throw cet::exception("RECO")<<"mu2e::RobustHelixFinder: stereo hit index not consistent"<< endl;

	  // see if this stereo hit has already been seen: if so, update its helix hit indices, if not create a map entry
	  auto ifnd = sthmap.find(stindex);
	  if (ifnd == sthmap.end()) {
	    std::pair<int, int> mypair = std::make_pair(-1,-1);
	    if (isfirst) mypair.first = ihh;
	    else         mypair.second = ihh;	         
	    sthmap[stindex] = mypair;
	  } else {
	    auto& mypair = ifnd->second;
	    if (isfirst) mypair.first = ihh;
	    else         mypair.second = ihh;
	  }
	}
      } 

      // now loop over the stereo hits in the map and update the positions
      for(auto imap=sthmap.begin();imap != sthmap.end(); ++imap) {
	const auto& sthit = sthcol.at(imap->first);

	CLHEP::Hep3Vector hdir;
	hseed.helix().direction(sthit.pos().z(),hdir);

	CLHEP::Hep3Vector pos1, pos2;
	sthit.position(shcol,tracker,pos1,pos2,hdir);

	if(imap->second.first >= 0) {
	  if (!retval) retval = hseed._hhits.at(imap->second.first).pos().diff2(pos1) > _dhit2;
	  hseed._hhits.at(imap->second.first)._pos = pos1;
	}
	if(imap->second.second >= 0) {
	  if (!retval)retval = hseed._hhits.at(imap->second.second).pos().diff2(pos2) > _dhit2;
	  hseed._hhits.at(imap->second.second)._pos = pos2;
	}
      }
    }
    return retval;
  }

}
using mu2e::RobustHelixFinder;
DEFINE_ART_MODULE(RobustHelixFinder);
