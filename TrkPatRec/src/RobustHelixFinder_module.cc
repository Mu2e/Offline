//
// Tracker Pattern Recognition based on Robust Helix Fit
//
// $Id: RobustHelixFinder_module.cc,v 1.2 2014/08/30 12:19:38 tassiell Exp $
// $Author: tassiell $
// $Date: 2014/08/30 12:19:38 $
//
// Original author D. Brown and G. Tassielli
//

#include "art/Framework/Principal/Event.h"
#include "GeometryService/inc/GeomHandle.hh"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
#include "GeneralUtilities/inc/Angles.hh"
#include "Mu2eUtilities/inc/MVATools.hh"

#include "DataProducts/inc/Helicity.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/StereoHit.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
#include "RecoDataProducts/inc/HelixSeed.hh"
#include "RecoDataProducts/inc/TrkFitFlag.hh"

#include "TrkReco/inc/TrkTimeCalculator.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"

#include "BTrk/BaBar/BaBar.hh"
#include "TrkReco/inc/RobustHelixFit.hh"
#include "TrkReco/inc/Chi2HelixFit.hh"

#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Mu2eUtilities/inc/ModuleHistToolBase.hh"
#include "Mu2eUtilities/inc/polyAtan2.hh"
#include "Mu2eUtilities/inc/HelixTool.hh"
#include "art/Utilities/make_tool.h"

#include "TrkPatRec/inc/RobustHelixFinder_types.hh"
#include "TrkReco/inc/RobustHelixFinderData.hh"
#include "TrkReco/inc/TrkFaceData.hh"

#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"

#include <boost/accumulators/accumulators.hpp>
#include "boost_fix/accumulators/statistics/stats.hpp"
#include "boost_fix/accumulators/statistics.hpp"
#include <boost/accumulators/statistics/median.hpp>

#include "TH1F.h"
#include "Math/VectorUtil.h"
#include "TVector2.h"

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
using namespace ROOT::Math::VectorUtil;

namespace {
  // comparison functor for sorting by z
  struct zcomp : public std::binary_function<mu2e::ComboHit,mu2e::ComboHit,bool> {
    bool operator()(mu2e::ComboHit const& p1, mu2e::ComboHit const& p2) { return p1._pos.z() < p2._pos.z(); }
  };

  // comparison functor for sorting byuniquePanel ID
  struct panelcomp : public std::binary_function<mu2e::ComboHit,mu2e::ComboHit,bool> {
    bool operator()(mu2e::ComboHit const& p1, mu2e::ComboHit const& p2) { return p1.strawId().uniquePanel() < p2.strawId().uniquePanel(); }
  };
  struct HelixHitMVA
  {
    std::vector <float> _pars,_pars2;
    float& _dtrans; // distance from hit to helix perp to the wrire
    float& _dwire;  // distance from hit to helix along the wrire
    float& _chisq;  // chisq of spatial information, using average errors
    float& _dt;     // time difference of hit WRT average
    float& _drho;   // hit transverse radius minus helix radius
    float& _dphi;   // hit azimuth minus helix azimuth (at the hit z)
    float& _rwdot;  // dot product between circle radial direction and wire direction
    float& _hrho;   // helix transverse radius (at the hit z)
    float& _hhrho;  // hit transverse radius
    //HelixHitMVA() : _pars(9,0.0),_dtrans(_pars[0]),_dwire(_pars[1]),_chisq(_pars[2]),_dt(_pars[3]),
    // _drho(_pars[4]),_dphi(_pars[5]),_rwdot(_pars[6]),_hrho(_pars[7]),_hhrho(_pars[8]){}
    HelixHitMVA() : _pars(7,0.0),_pars2(2,0.0),_dtrans(_pars[0]),_dwire(_pars[1]),_chisq(_pars[2]),_dt(_pars[3]),
		    _drho(_pars[4]),_dphi(_pars[5]),_rwdot(_pars[6]),_hrho(_pars[0]),_hhrho(_pars2[1]) {}
  };

}

namespace mu2e {

  class RobustHelixFinder : public art::EDProducer {
  public:
    explicit RobustHelixFinder(fhicl::ParameterSet const&);
    virtual ~RobustHelixFinder();
    virtual void beginJob();
    virtual void beginRun(art::Run&   run   );
    virtual void produce(art::Event& event );

  private:
    int                                 _diag,_debug,_reducedchi2;
    int                                 _printfreq;
    bool				_prefilter; // prefilter hits based on sector
    bool				_updatestereo; // update the stereo hit positions each iteration
    int 				_minnsh; // minimum # of strawHits to work with
    float			        _pitch; // average pitch to assume in time calculations
    float                               _maxchi2dxy;
    float                               _maxchi2dzphi;
    float                               _maxphihitchi2;
    float				_maxdr; // maximum hit-helix radius difference
    float				_maxrpull; // maximum hit-helix radius difference pull
    bool                                _targetconInit;//require the firs circle fit to intersect the Al stopping Target
    bool                                _targetcon;//require the circle fit to intersect the Al stopping Target
    float                               _rpullScaleF;//need to scale the radial pull in filterCircleHits
    float				_maxphisep; // maximum separation in global azimuth of hits
    TrkFitFlag				_saveflag; // write out all helices that satisfy these flags
    unsigned				_maxniter;  // maximum # of iterations over outlier filtering + fitting
    float				_cradres; // average center resolution along center position (mm)
    float				_cperpres; // average center resolution perp to center position (mm)
    float				_maxdwire; // outlier cut on distance between hit and helix along wire
    float				_maxdtrans; // outlier cut on distance between hit and helix perp to wire
    float				_maxchisq; // outlier cut on chisquared
    float				_maxrwdot; // outlier cut on angle between radial direction and wire: smaller is better
    float				_minrerr; // minimum radius error

    bool				_usemva; // use MVA to cut outliers
    float                               _minmva; // outlier cut on MVA

    art::ProductToken<ComboHitCollection> const _chToken;
    art::ProductToken<TimeClusterCollection> const _tcToken;

    StrawHitFlag  _hsel, _hbkg;

    MVATools _stmva, _nsmva;
    HelixHitMVA _vmva; // input variables to TMVA for filtering hits

    TH1F* _niter, *_niterxy, *_niterfz, *_nitermva;

    RobustHelixFit   _hfit;
    Chi2HelixFit     _chi2hfit;

    std::vector<Helicity> _hels; // helicity values to fit
    TrkTimeCalculator _ttcalc;
    StrawHitFlag      _outlier;
    bool              _updateStereo;
    
    std::unique_ptr<ModuleHistToolBase>   _hmanager;
    RobustHelixFinderTypes::Data_t        _data;
    RobustHelixFinderData                 _hfResult;

    void     findHelices(ComboHitCollection& chcol, const TimeClusterCollection& tccol);    
    void     prefilterHits(RobustHelixFinderData& helixData, int& nFilteredStrawHits); 
    unsigned filterCircleHits(RobustHelixFinderData& helixData); 
    int      filterChi2ZPhiHits(RobustHelixFinderData& helixData);
    int      filterChi2XYHits(RobustHelixFinderData& helixData); 
    bool     filterHits(RobustHelixFinderData& helixData);
    void     fillMVA(RobustHelixFinderData& helixData); 
    bool     filterHitsMVA(RobustHelixFinderData& helixData);
    void     updateT0(RobustHelixFinderData& helixData);
    bool     updateStereo(RobustHelixFinderData& helixData);
    unsigned hitCount(RobustHelixFinderData& helixData);
    void     pickBestHelix      (std::vector<HelixSeed>& HelVec, int &Index_best);
    void     fillFaceOrderedHits(RobustHelixFinderData& helixData);
    void     fillGoodHits       (RobustHelixFinderData& helixData);
    void     fitHelix           (RobustHelixFinderData& helixData);
    void     fitChi2Helix       (RobustHelixFinderData& helixData);
    void     refitHelix         (RobustHelixFinderData& helixData);
    void     findMissingHits    (RobustHelixFinderData& helixData);
    void     fillPluginDiag     (RobustHelixFinderData& helixData, int helCounter);
    void     updateChi2HelixInfo(RobustHelixFinderData& helixData);
    void     updateHelixXYInfo  (RobustHelixFinderData& helixData);
    void     updateHelixZPhiInfo(RobustHelixFinderData& helixData);
    void     searchWorstHitXY   (RobustHelixFinderData& helixData, HitInfo_t& hitInfo);
    void     searchWorstHitZPhi (RobustHelixFinderData& helixData, HitInfo_t& hitInfo);
  };

  RobustHelixFinder::RobustHelixFinder(fhicl::ParameterSet const& pset) :
    art::EDProducer{pset},
    _diag        (pset.get<int>("diagLevel",0)),
    _debug       (pset.get<int>("debugLevel",0)),
    _reducedchi2 (pset.get<int>("reducedchi2",0)),
    _printfreq   (pset.get<int>("printFrequency",101)),
    _prefilter   (pset.get<bool>("PrefilterHits",true)),
    _updatestereo(pset.get<bool>("UpdateStereoHits",false)),
    _minnsh      (pset.get<int>("minNStrawHits",10)),
    _pitch       (pset.get<int>("AveragePitch",0.6)),
    _maxchi2dxy  (pset.get<float>("MaxChi2dXY", 5.0)),
    _maxchi2dzphi(pset.get<float>("MaxChi2dZPhi", 5.0)),
    _maxphihitchi2(pset.get<float>("MaxHitPhiChi2", 25.0)),
    _maxdr	 (pset.get<float>("MaxRadiusDiff",100.0)), // mm
    _maxrpull	 (pset.get<float>("MaxRPull",5.0)), // unitless
    _targetconInit(pset.get<bool>("targetconsistent_init",true)),
    _targetcon   (pset.get<bool>("targetconsistent",true)),
    _rpullScaleF (pset.get<float>("RPullScaleF",1.414)), // unitless
    _maxphisep	 (pset.get<float>("MaxPhiHitSeparation",1.0)),
    _saveflag    (pset.get<vector<string> >("SaveHelixFlag",vector<string>{"HelixOK"})),
    _maxniter    (pset.get<unsigned>("MaxIterations",10)), // iterations over outlier removal
    _cradres     (pset.get<float>("CenterRadialResolution",20.0)),
    _cperpres    (pset.get<float>("CenterPerpResolution",12.0)),
    _maxdwire    (pset.get<float>("MaxWireDistance",200.0)), // max distance along wire
    _maxdtrans   (pset.get<float>("MaxTransDistance",80.0)), // max distance perp to wire (and z)
    _maxchisq    (pset.get<float>("MaxChisquared", 25.)), //100.0)), // max chisquared
    _maxrwdot	 (pset.get<float>("MaxRWDot",1.0)),
    _minrerr     (pset.get<float>("MinRadiusErr",20.0)), // mm
    _usemva      (pset.get<bool>("UseHitMVA",false)),
    _minmva      (pset.get<float> ("MinMVA",0.1)), // min MVA output to define an outlier
    _chToken{consumes<ComboHitCollection>(pset.get<art::InputTag>("ComboHitCollection"))},
    _tcToken{consumes<TimeClusterCollection>(pset.get<art::InputTag>("TimeClusterCollection"))},
    _hsel        (pset.get<std::vector<std::string> >("HitSelectionBits",std::vector<string>{"TimeDivision"})),
    _hbkg        (pset.get<std::vector<std::string> >("HitBackgroundBits",std::vector<std::string>{"Background"})),
    _stmva       (pset.get<fhicl::ParameterSet>("HelixStereoHitMVA",fhicl::ParameterSet())),
    _nsmva       (pset.get<fhicl::ParameterSet>("HelixNonStereoHitMVA",fhicl::ParameterSet())),
    _hfit        (pset.get<fhicl::ParameterSet>("RobustHelixFit",fhicl::ParameterSet())),
    _chi2hfit    (pset.get<fhicl::ParameterSet>("Chi2HelixFit",fhicl::ParameterSet())),
    _ttcalc      (pset.get<fhicl::ParameterSet>("T0Calculator",fhicl::ParameterSet())),
    _outlier     (StrawHitFlag::outlier),
    _updateStereo    (pset.get<bool>("UpdateStereo",false))
  {
    std::vector<int> helvals = pset.get<std::vector<int> >("Helicities",vector<int>{Helicity::neghel,Helicity::poshel});
    for(auto hv : helvals) {
      Helicity hel(hv);
      _hels.push_back(hel);
      produces<HelixSeedCollection>(Helicity::name(hel));
    }
    
     if (_diag != 0) _hmanager = art::make_tool<ModuleHistToolBase>(pset.get<fhicl::ParameterSet>("diagPlugin"));
     else            _hmanager = std::make_unique<ModuleHistToolBase>();

     //    _data.result    = &_hfit;

    _chi2hfit.setRobustHelixFitter(&_hfit);

  }
  
  RobustHelixFinder::~RobustHelixFinder(){}

  //-----------------------------------------------------------------------------
  void RobustHelixFinder::beginRun(art::Run& ) {
    mu2e::GeomHandle<mu2e::Tracker> th;
    const Tracker* tracker = th.get();

    mu2e::GeomHandle<mu2e::Calorimeter> ch;

    _hfit.setTracker    (tracker);
    _hfit.setCalorimeter(ch.get());
    
    _chi2hfit.setTracker    (tracker);
    _chi2hfit.setCalorimeter(ch.get());
  }
  //--------------------------------------------------------------------------------

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
      _niterxy = tfs->make<TH1F>( "niterxy" , "Number of XY Fit Iteraions",201,-0.5,200.5);
      _niterfz = tfs->make<TH1F>( "niterfz" , "Number of FZ Fit Iteraions",201,-0.5,200.5);
      _nitermva = tfs->make<TH1F>( "nitermva" , "Number of MVA Fit Iteraions",201,-0.5,200.5);
      _hmanager->bookHistograms(tfs);
    }
  }

  void RobustHelixFinder::produce(art::Event& event ) {
    // find input
    auto const& tcH = event.getValidHandle(_tcToken);
    const TimeClusterCollection& tccol(*tcH);

    auto const& chH = event.getValidHandle(_chToken);
    const ComboHitCollection& chcol(*chH);

    // create output: seperate by helicity
    std::map<Helicity,unique_ptr<HelixSeedCollection>> helcols;
    int counter(0);
    for( auto const& hel : _hels) {
      helcols[hel] = unique_ptr<HelixSeedCollection>(new HelixSeedCollection());
      _data.nseeds [counter] = 0;
      ++counter;
    }
    
    _data.event       = &event;
    //    _data.result      = &_hfit;
    _data.nTimePeaks  = tccol.size();

    _hfResult._chcol  = &chcol;
      
    // create initial helicies from time clusters: to begin, don't specificy helicity
    for (size_t index=0;index< tccol.size();++index) {
      const auto& tclust = tccol[index];
      HelixSeed hseed;
      hseed._status.merge(TrkFitFlag::TPRHelix);
      //clear the variables in hfResult
      _hfResult.clearTempVariables();

      //set variables used for searching the helix candidate
      _hfResult._hseed              = hseed;
      _hfResult._timeCluster        = &tclust;
      _hfResult._hseed._hhits.setParent(chcol.parent());
      _hfResult._hseed._t0          = tclust._t0;
      _hfResult._hseed._timeCluster = art::Ptr<TimeCluster>(tcH,index);
      // copy combo hits
      fillFaceOrderedHits(_hfResult);

      //skip the reconstruction if there are few strawHits
      if (_hfResult._nFiltStrawHits < _minnsh)                  continue;

      // filter hits and test
      int nFilteredSh(0);
      if (_prefilter) prefilterHits(_hfResult,nFilteredSh);

      if ((_hfResult._nFiltStrawHits - nFilteredSh) < _minnsh)  continue;   

      _hfResult._hseed._status.merge(TrkFitFlag::hitsOK);
      if (_diag) _hfResult._diag.circleFitCounter = 0;

      // initial circle fit

      if (_reducedchi2){
	_chi2hfit.fitChi2Circle(_hfResult, _targetcon);
      }else{
	_hfit.fitCircle(_hfResult, _targetconInit);//require consistency for the trajectory of being produced in the Al stopping target
      }

      if (_diag && _reducedchi2) {
	_hfResult._diag.nShFitCircle = _hfResult._nXYSh;
	_hfResult._diag.nChFitCircle = _hfResult._sxy.qn()-1;//take into account one hit form the stopping target center
      }
      //check the number of points associated with the result of the circle fit
      // if (_hfResult._nXYSh < _minnsh)                           continue;
      
      if (_hfResult._hseed._status.hasAnyProperty(TrkFitFlag::circleOK)) {
	// loop over helicities. 
	unsigned    helCounter(0);
	HelixSeed   helixSeed_from_fitCircle = _hfResult._hseed;

	std::vector<HelixSeed>          helix_seed_vec;
      
	for(auto const& hel : _hels ) {
	  // tentatively put a copy with the specified helicity in the appropriate output vector
	  RobustHelixFinderData tmpResult(_hfResult);
	  tmpResult._hseed._helix._helicity = hel;

	  //fit the helix: refine the XY-circle fit + performs the ZPhi fit
	  // it also performs a clean-up of the hits with large residuals
	  if (_reducedchi2)
	    fitChi2Helix(tmpResult);
	  else
	    fitHelix(tmpResult);
	    

	  if (tmpResult._hseed.status().hasAnyProperty(_saveflag)){
	    //fill the hits in the HelixSeedCollection
	    fillGoodHits(tmpResult);
	    
	    helix_seed_vec.push_back(tmpResult._hseed);

	    // HelixSeedCollection* hcol = helcols[hel].get();
	    // hcol->push_back(tmpResult._hseed);

	    if (_diag > 0) {
	      fillPluginDiag(tmpResult, helCounter);
	    }
	  }
	  ++helCounter;
	}//end loop over the helicity

	if (helix_seed_vec.size() == 0)                       continue;

	int    index_best(-1);
	pickBestHelix(helix_seed_vec, index_best);

	if ( (index_best>=0) && (index_best < 2) ){
	  Helicity              hel_best = helix_seed_vec[index_best]._helix._helicity;
	  HelixSeedCollection*  hcol     = helcols[hel_best].get();
	  hcol->push_back(helix_seed_vec[index_best]);
	} else if (index_best == 2){//both helices need to be saved
	
	  for (unsigned k=0; k<_hels.size(); ++k){
	    Helicity              hel   = helix_seed_vec[k]._helix._helicity;
	    HelixSeedCollection*  hcol  = helcols[hel].get();
	    hcol->push_back(helix_seed_vec[k]);
	  }
	}	

      }	
      
    }
    // put final collections into event 
    if (_diag > 0) _hmanager->fillHistograms(&_data);

    for(auto const& hel : _hels ) {
      event.put(std::move(helcols[hel]),Helicity::name(hel));
    }
  }
//--------------------------------------------------------------------------------
// function to select the best Helix among the results of the two helicity hypo
//--------------------------------------------------------------------------------
  void  RobustHelixFinder::pickBestHelix(std::vector<HelixSeed>& HelVec, int &Index_best){
    if (HelVec.size() == 1) {
      Index_best = 0;
      return;
    }
    
    const HelixSeed           *h1, *h2;
    const ComboHitCollection  *tlist, *clist;
    int                        nh1, nh2;

    h1     = &HelVec[0];
//------------------------------------------------------------------------------
// check if an AlgorithmID collection has been created by the process
//-----------------------------------------------------------------------------
    tlist  = &h1->hits();
    nh1    = tlist->size();

    h2     = &HelVec[1];
//-----------------------------------------------------------------------------
// at Mu2e, 2 helices with different helicity could be duplicates of each other
//-----------------------------------------------------------------------------
    clist  = &h2->hits();
    nh2    = clist->size();
//-----------------------------------------------------------------------------
// pick the helix with the largest number of hits
//-----------------------------------------------------------------------------
    if (nh2 > nh1) {
//-----------------------------------------------------------------------------
// h2 is a winner, no need to save h1
//-----------------------------------------------------------------------------
      Index_best = 1;
      return;
    }
    else if (nh1 > nh2){
//-----------------------------------------------------------------------------
// h1 is a winner, mark h2 in hope that it will be OK, continue looping
//-----------------------------------------------------------------------------
      Index_best = 0;
      return;
    }
    
//-----------------------------------------------------------------------------
// in case they have the exact amount of hits, pick the one with better chi2dZphi
//-----------------------------------------------------------------------------
    if (nh1 == nh2) {
      float   chi2dZphi_h1 = h1->helix().chi2dZPhi();
      float   chi2dZphi_h2 = h2->helix().chi2dZPhi();
      if (chi2dZphi_h1 < chi2dZphi_h2){
	Index_best = 0;
	return;
      }else {
	Index_best = 1;
	return;      
      }
    }

  }

//--------------------------------------------------------------------------------
// 
//--------------------------------------------------------------------------------
  void RobustHelixFinder::fillGoodHits(RobustHelixFinderData& helixData){
    
    ComboHit*     hit(0);
    unsigned      nhits = helixData._chHitsToProcess.size();

    for (unsigned f=0; f<nhits; ++f){
      hit = &helixData._chHitsToProcess[f];
      if (hit->_flag.hasAnyProperty(_outlier))     continue;
      
      ComboHit                hhit(*hit);					
      helixData._hseed._hhits.push_back(hhit);
    }

    if (_diag){
      HelixTool helTool(&helixData._hseed, 3);
      helixData._diag.nLoops            = helTool.nLoops();
      helixData._diag.meanHitRadialDist = helTool.meanHitRadialDist();
    }
  }


  void RobustHelixFinder::fillMVA(RobustHelixFinderData& helixData)
  {
    RobustHelix& helix = helixData._hseed._helix;

    static XYZVec  zaxis(0.0,0.0,1.0); // unit in z direction
    ComboHit*      hhit(0);
    
    for (unsigned f=0; f<helixData._chHitsToProcess.size(); ++f){
      hhit = &helixData._chHitsToProcess[f];
	
      if (hhit->_flag.hasAnyProperty(_outlier))   continue;

      const XYZVec& wdir = hhit->wdir();
      XYZVec wtdir = zaxis.Cross(wdir); // transverse direction to the wire
      XYZVec cvec = PerpVector(hhit->pos() - helix.center(),Geom::ZDir());// direction from the circle center to the hit
      XYZVec cdir = cvec.Unit();        // direction from the circle center to the hit
      XYZVec cperp = zaxis.Cross(cdir); // direction perp to the radius

      XYZVec hpos = hhit->pos();      // this sets the z position to the hit z
      helix.position(hpos);                     // this computes the helix expectation at that z
      XYZVec dh = hhit->pos() - hpos; // this is the vector between them

      _vmva._dtrans = fabs(dh.Dot(wtdir));              // transverse projection
      _vmva._dwire = fabs(dh.Dot(wdir));               // projection along wire direction
      _vmva._drho = fabs(sqrtf(cvec.mag2()) - helix.radius()); // radius difference
      _vmva._dphi = fabs(hhit->helixPhi() - helix.circleAzimuth(hhit->pos().z())); // azimuth difference WRT circle center
      _vmva._hhrho = sqrtf(cvec.mag2());            // hit transverse radius WRT circle center
      _vmva._hrho = sqrtf(hpos.Perp2());            // hit detector transverse radius
      _vmva._rwdot = fabs(wdir.Dot(cdir));  // compare directions of radius and wire

      // compute the total resolution including hit and helix parameters first along the wire
      float wres2 = std::pow(hhit->posRes(StrawHitPosition::wire),(int)2) +
	std::pow(_cradres*cdir.Dot(wdir),(int)2) +
	std::pow(_cperpres*cperp.Dot(wdir),(int)2);

      // transverse to the wires
      float wtres2 = std::pow(hhit->posRes(StrawHitPosition::trans),(int)2) +
	std::pow(_cradres*cdir.Dot(wtdir),(int)2) +
	std::pow(_cperpres*cperp.Dot(wtdir),(int)2);

      _vmva._chisq = sqrtf( _vmva._dwire*_vmva._dwire/wres2 + _vmva._dtrans*_vmva._dtrans/wtres2 );          
      _vmva._dt = hhit->time() - helixData._hseed._t0.t0();

      if (hhit->_flag.hasAnyProperty(StrawHitFlag::stereo))
	{
	  hhit->_qual = _stmva.evalMVA(_vmva._pars);
	} else {
	hhit->_qual = _nsmva.evalMVA(_vmva._pars);
      }
    }
  }

  bool RobustHelixFinder::filterHitsMVA(RobustHelixFinderData& helixData)
  {  
    bool           changed(false);
    ComboHit*      hhit(0);
    
    for (unsigned f=0; f<helixData._chHitsToProcess.size(); ++f){
 
      hhit = &helixData._chHitsToProcess[f];
	
      bool oldout = hhit->_flag.hasAnyProperty(_outlier);

      if (hhit->_qual < _minmva ) hhit->_flag.merge(_outlier);
      else                        hhit->_flag.clear(_outlier);

      changed |= oldout != hhit->_flag.hasAnyProperty(_outlier);
    }    

    return changed;
  }

  int  RobustHelixFinder::filterChi2XYHits(RobustHelixFinderData& helixData)
  {
    //reset the value of the XY fit result
    helixData._hseed._status.clear(TrkFitFlag::circleOK);

  
    ComboHit*     hit(0);

    //perform a reduced chi2 fit
    _chi2hfit.refineFitXY(helixData, _targetcon);

    int           changed(0);
    int           oldNHitsSh = helixData._nXYSh;
    
    RobustHelix&  helix      = helixData._hseed._helix;

    // helixData._hseed._status.clear(TrkFitFlag::circleOK);

    if (helixData._nXYSh >= _minnsh) {//update the helix info
      //need to update the weights in the LSqsum
      _chi2hfit.refineFitXY(helixData, _targetcon);
   
      //      updateHelixXYInfo(helixData);//should be unnecessary!FIXME!
          
      //search and remove the worst hit(s) if necessary
      HitInfo_t  worstHit;
      float      chi2d = helixData._sxy.chi2DofCircle();
      
      while( (chi2d > _maxchi2dxy) && (helixData._nXYSh >= _minnsh) && (worstHit.face >=0)){
	//reset the content of the worstHit
	worstHit.face          = -1;
	worstHit.panel         = -1;
	worstHit.panelHitIndex = -1;
	//	worstHit.weightXY      =  0;

	searchWorstHitXY(helixData, worstHit);
	
	if (worstHit.face >=0){//check if a bad was found or not
	  hit    = &helixData._chHitsToProcess[worstHit.panelHitIndex];
	  
	  hit->_flag.merge(_outlier);

	  helixData._sxy.removePoint(hit->pos().x(), hit->pos().y(), hit->_xyWeight);//worstHit.weightXY);
	  helixData._nXYSh -= hit->nStrawHits();
	  helixData._nXYCh -= 1;
	  _chi2hfit.refineFitXY(helixData, _targetcon);//should be unnecessary!FIXME!
	  chi2d             = helixData._sxy.chi2DofCircle();
	}
      }

      //at this point
      if (_hfit.goodCircle(helix) && (helixData._nXYSh >= _minnsh))  {
	helixData._hseed._status.merge(TrkFitFlag::circleOK);

	if (_diag){
	  helixData._diag.rsxy_1     = helix._radius;
	  helixData._diag.chi2dsxy_1 = helixData._sxy.chi2DofCircle();
	  helixData._diag.nshsxy_1   = helixData._nXYSh;
	}
      }
    }
 
    changed = oldNHitsSh - helixData._nXYSh;

    return changed;
  }


  // 3d selection on top of radial selection
  int RobustHelixFinder::filterChi2ZPhiHits(RobustHelixFinderData& helixData)
  {

    //check if the initial value of lambda and phi0 are physical
    if (!helixData._hseed._status.hasAnyProperty(TrkFitFlag::phizOK))  return false;

    //reset the value of the ZPhi fit result
    helixData._hseed._status.clear(TrkFitFlag::phizOK);
    
    RobustHelix&  helix  = helixData._hseed._helix;

    int           changed(0);
    int           oldNHitsSh = helixData._nZPhiSh;

    ComboHit*     hit(0);

    // float         z, phi, phi_ref, dx, dy, dphi, resid, wt;

    //    helixData._hseed._status.clear(TrkFitFlag::circleOK);
    _chi2hfit.refineFitZPhi(helixData);

    if (helixData._nZPhiSh >= _minnsh) {//update the helix info
      //need to update the weights in the LSqsum
      _chi2hfit.refineFitZPhi(helixData);

      //      updateHelixZPhiInfo(helixData);//should be unnecessary!FIXME!

      //search and remove the worst hit(s) if necessary
      HitInfo_t  worstHit;
      float      chi2d = helixData._szphi.chi2DofLine();
      
      while( (chi2d > _maxchi2dzphi) && (helixData._nZPhiSh >= _minnsh) && (worstHit.face >=0)){
	//reset the content of the worstHit
	worstHit.face          = -1;
	worstHit.panel         = -1;
	worstHit.panelHitIndex = -1;
	
	searchWorstHitZPhi(helixData, worstHit);
	
	if (worstHit.face >=0){//check if a bad was found or not
	  hit    = &helixData._chHitsToProcess[worstHit.panelHitIndex];
	  
	  hit->_flag.merge(_outlier);

	  helixData._szphi.removePoint(hit->pos().z(), hit->helixPhi(), hit->_zphiWeight);
	  helixData._nZPhiSh -= hit->nStrawHits();
	  _chi2hfit.refineFitZPhi(helixData);
	  //	  updateHelixZPhiInfo(helixData);//should be unnecessary!FIXME!
	  chi2d               = helixData._szphi.chi2DofLine();
	}
      }

      //at this point
      if (_hfit.goodFZ(helix) && (helixData._nZPhiSh >= _minnsh))  {
	helixData._hseed._status.merge(TrkFitFlag::phizOK);

	if (_diag){
	  helixData._diag.lambdaszphi_1 = helix._lambda;
	  helixData._diag.chi2dszphi_1  = helixData._szphi.chi2DofLine();
	  helixData._diag.nshszphi_1    = helixData._nZPhiSh;	
	}
      }
    }

    changed = oldNHitsSh - helixData._nZPhiSh;

    return changed;
  }


  // 3d selection on top of radial selection
  bool RobustHelixFinder::filterHits(RobustHelixFinderData& helixData)
  {
    RobustHelix& helix = helixData._hseed._helix;
    bool changed(false);
    static XYZVec zaxis(0.0,0.0,1.0); // unit in z direction
    int      nGoodSH(0);

    // loop over hits
    ComboHit*     hit(0);
    FaceZ_t*      facez;

    int           nhitsFace(0);
    float         chCounter(1e-10), chi2dZPhi(0);

    for (int f=0; f<StrawId::_ntotalfaces; ++f){
      facez     = &helixData._oTracker[f];
      
      float      minChi2(_maxchisq);
      HitInfo_t  indexBestComboHit;

      nhitsFace = facez->nChHits();
      if (nhitsFace == 0)                        continue;
      int        idFirstFaceCh(facez->idChBegin);
      for (int ip=0; ip<nhitsFace; ++ip){
	hit = &helixData._chHitsToProcess[idFirstFaceCh + ip];
	bool trash=hit->_flag.hasAnyProperty(_outlier);
	if (trash)                               continue;

	float hphi = polyAtan2(hit->pos().y(),hit->pos().x());//phi();
	float dphi = fabs(Angles::deltaPhi(hphi,helix.fcent()));

	const XYZVec& wdir = hit->wdir();
	XYZVec wtdir = zaxis.Cross(wdir);   // transverse direction to the wire
	XYZVec cvec = PerpVector(hit->pos() - helix.center(),Geom::ZDir()); // direction from the circle center to the hit
	XYZVec cdir = cvec.Unit();          // direction from the circle center to the hit
	XYZVec cperp = zaxis.Cross(cdir);   // direction perp to the radius

	XYZVec hpos = hit->pos(); // this sets the z position to the hit z
	helix.position(hpos);                // this computes the helix expectation at that z
	XYZVec dh = hit->pos() - hpos;   // this is the vector between them
	float dtrans = fabs(dh.Dot(wtdir)); // transverse projection
	float dwire = fabs(dh.Dot(wdir));   // projection along wire direction

	// compute the total resolution including hit and helix parameters first along the wire
	float wres2 = std::pow(hit->posRes(StrawHitPosition::wire),(int)2) +
	  std::pow(_cradres*cdir.Dot(wdir),(int)2) +
	  std::pow(_cperpres*cperp.Dot(wdir),(int)2);
	// transverse to the wires
	float wtres2 = std::pow(hit->posRes(StrawHitPosition::trans),(int)2) +
	  std::pow(_cradres*cdir.Dot(wtdir),(int)2) +
	  std::pow(_cperpres*cperp.Dot(wtdir),(int)2);

	float chisq = dwire*dwire/wres2 + dtrans*dtrans/wtres2;
	
	if( dphi > _maxphisep || fabs(dwire) > _maxdwire || fabs(dtrans) > _maxdtrans || chisq > _maxchisq) 
	  {
	    changed = true;
	  }
 
	if ( chisq <= minChi2)
	  {
	    minChi2 = chisq;
	      
	    indexBestComboHit.face          = f;
	    indexBestComboHit.panel         = hit->strawId().uniquePanel();
	    indexBestComboHit.panelHitIndex = facez->idChBegin + ip;
	  }

	//flagg all hits within the face as outlier. Only the best found will be "cleared"
	hit->_flag.merge(_outlier);

      }//end loop over the panels
      
      //remove the outlier flag 
      if (indexBestComboHit.face >=0 ) {
	hit     = &helixData._chHitsToProcess[indexBestComboHit.panelHitIndex];

	//remove the outlier flag
	hit->_flag.clear(StrawHitFlag::outlier);
	nGoodSH += hit->nStrawHits();
	
	chi2dZPhi += minChi2;
	chCounter += 1.;
      }
    }//end loop over the faces
    
    helixData._nZPhiSh = nGoodSH;

    //update the value of the chi2ZPhi
    helix._chi2dZPhi   = chi2dZPhi/chCounter;

    if (_diag) {
      helixData._diag.chi2dZPhi = chi2dZPhi/chCounter;
    }

    return changed;
  }
    
  void     RobustHelixFinder::findMissingHits(RobustHelixFinderData& helixData){
    FaceZ_t*   facez;

    ComboHit*  hit(0);
    int        nhitsPerPanel(0), n_added_points(0);
    HitInfo_t  bestHit;
    
    float      wtXY(0),wtZPhi(0),dr(0),dphi(0), phi_pred(0);
    float      drChi2, dphiChi2, hitChi2Max(_maxchi2dxy), hitChi2;

    //get  dfdz and phi0
    float      dfdz = helixData._szphi.dfdz();
    float      phi0 = helixData._szphi.phi0();
    //get the circle info
    float      r    = helixData._sxy.radius();
    XYVec      helCenter;
    helCenter.SetX( helixData._sxy.x0());
    helCenter.SetY( helixData._sxy.y0());


  NEXT_ITERATION:;
    //reset the info of the best-hit
    bestHit.face          = -1;	 
    bestHit.panel         = -1;
    bestHit.panelHitIndex = -1;
    // bestHit.weightXY      = 1.;
    // bestHit.weightZPhi    = 1.;
    
    hitChi2Max            = _maxchi2dxy;
    
    for (int f=0; f<StrawId::_ntotalfaces; ++f){
      facez         = &helixData._oTracker[f];
      bool       isFaceUsed(false);
      HitInfo_t  hitInfo;

      nhitsPerPanel  = facez->nChHits();
	  
      if (nhitsPerPanel == 0)                       continue;
      
      for (int i=0; i<nhitsPerPanel;++i){
	hit =  &helixData._chHitsToProcess[facez->idChBegin +i];
	if (!hit->_flag.hasAnyProperty(_outlier))   {
	  isFaceUsed = true;
	  break;//skip the faces where there is already a hit
	}
	XYVec rvec = (XYVec(hit->pos().x(),hit->pos().y())-helCenter);
	dr       = sqrtf(rvec.Mag2()) - r;
	wtXY     = _chi2hfit.evalWeightXY(*hit, helCenter);
	drChi2   = sqrtf(dr*dr*wtXY);

	phi_pred = hit->pos().z()*dfdz + phi0;
	dphi     = phi_pred - hit->helixPhi();
	wtZPhi   = _chi2hfit.evalWeightZPhi(*hit,helCenter,r);
	dphiChi2 = sqrtf(dphi*dphi*wtZPhi);
	
	hitChi2  = (drChi2 + dphiChi2)/2.;
	
	if ( (drChi2<_maxchi2dxy) && (dphiChi2<_maxchi2dzphi) && (hitChi2 < hitChi2Max)){
	  hitChi2Max  = hitChi2;
	  
	  hitInfo.face          = f;	 
	  hitInfo.panel         = hit->strawId().uniquePanel();	 
	  hitInfo.panelHitIndex = facez->idChBegin + i;

	  //update weight info
	  hit->_xyWeight        = wtXY;
	  hit->_zphiWeight      = wtZPhi;
   
	}
      }//end loop over the hits within the panel

      if( (!isFaceUsed) && (hitInfo.face >=0 ) ) {
	bestHit.face          = hitInfo.face         ;	 
	bestHit.panel         = hitInfo.panel        ;	 
	bestHit.panelHitIndex = hitInfo.panelHitIndex;
      }
    }//end loop pver the faces
    
    
    if ( (bestHit.face >= 0) ){
      hit    = &helixData._chHitsToProcess[bestHit.panelHitIndex];
	

      //add the point 
      hit->_flag.clear(_outlier);

      helixData._sxy.addPoint(hit->pos().x(), hit->pos().y(), hit->_xyWeight);//bestHit.weightXY);
      helixData._nXYSh   += hit->nStrawHits();

      helixData._szphi.addPoint(hit->pos().z(), hit->helixPhi(), hit->_zphiWeight);//bestHit.weightZPhi);
      helixData._nZPhiSh += hit->nStrawHits();


      //update the helix
      updateChi2HelixInfo(helixData);
      
      ++n_added_points;
	                              goto NEXT_ITERATION;
    }
      
    if (_diag) helixData._diag.nrescuedhits = n_added_points;

  }




  void RobustHelixFinder::prefilterHits(RobustHelixFinderData& HelixData, int& NRemovedStrawHits)
  {
    // ComboHitCollection& hhits = HelixData._hseed._hhits;

    bool changed(true);
    // size_t nhit = hhits.size();
    int nhit = HelixData._nFiltComboHits;

    ComboHit*  hit(0);
    ComboHit*  worsthit(0);

    NRemovedStrawHits = 0;

    while (changed && nhit > 0)
      {
	nhit = 0;
	changed = false;
	accumulator_set<float, stats<tag::median(with_p_square_quantile) > > accx;
	accumulator_set<float, stats<tag::median(with_p_square_quantile) > > accy;

	for (unsigned f=0; f<HelixData._chHitsToProcess.size(); ++f){
	  hit =  &HelixData._chHitsToProcess[f];
	  bool trashHit=hit->_flag.hasAnyProperty(_outlier);
	  if (trashHit)                              continue;
	  accx(hit->_pos.x());
	  accy(hit->_pos.y());
	  ++nhit;
	}
	
	float mx = extract_result<tag::median>(accx);
	float my = extract_result<tag::median>(accy);
	float mphi = polyAtan2(my,mx);//atan2f(my,mx);

	float maxdphi{0.0};
	// auto worsthit = hhits.end();
	for (unsigned f=0; f<HelixData._chHitsToProcess.size(); ++f){
	  hit =  &HelixData._chHitsToProcess[f];
	  bool trashHit = hit->_flag.hasAnyProperty(_outlier);
	  if (trashHit)                              continue;
	  float phi  = polyAtan2(hit->pos().y(), hit->pos().x());//ihit->pos().phi();
	  float dphi = fabs(Angles::deltaPhi(phi,mphi));
	  if(dphi > maxdphi)
	    {
	      maxdphi = dphi;
	      worsthit = hit;
	    }
	}//end loop over the faces

	if (maxdphi > _maxphisep)
	  {
	    worsthit->_flag.merge(_outlier);
	    NRemovedStrawHits += worsthit->nStrawHits();
	    changed = true;
	  }
      }
  }

  void RobustHelixFinder::updateT0(RobustHelixFinderData& helixData)
  {
  // compute the pitch
    float pitch = helixData._hseed.helix().pitch();
    accumulator_set<float, stats<tag::weighted_variance(lazy)>, float > terr;
  // update t0 from calo cluster according to current pitch 
    if (helixData._hseed.caloCluster().isNonnull()){
      float cwt = std::pow(1.0/_ttcalc.caloClusterTimeErr(),2);
      terr(_ttcalc.caloClusterTime(*helixData._hseed.caloCluster(),pitch),weight=cwt);
    }
    ComboHit*      hit(0);
    float hwt = std::pow(1.0/_ttcalc.strawHitTimeErr(),2);
    for (unsigned f=0; f<helixData._chHitsToProcess.size(); ++f){
      hit = &helixData._chHitsToProcess[f];
      if (hit->_flag.hasAnyProperty(_outlier))   continue;
      terr(_ttcalc.comboHitTime(*hit,pitch),weight=hwt);
    }//end faces loop
    if (sum_of_weights(terr) > 0.0)
      {
	helixData._hseed._t0._t0 = extract_result<tag::weighted_mean>(terr);
	helixData._hseed._t0._t0err = sqrtf(std::max(float(0.0),extract_result<tag::weighted_variance(lazy)>(terr))/extract_result<tag::count>(terr));
      }
  }

  //------------------------------------------------------------------------------------------
  void     RobustHelixFinder::fillFaceOrderedHits(RobustHelixFinderData& HelixData){
  
    const vector<StrawHitIndex>& shIndices = HelixData._timeCluster->hits();
    mu2e::RobustHelixFinderData::ChannelID cx, co;

    int     size  = shIndices.size();
    int     nFiltComboHits(0), nFiltStrawHits(0);
    // int nTotalStations = _tracker->nStations();
    //--------------------------------------------------------------------------------
    int loc;
    StrawHitFlag flag;

    //sort the hits by z coordinate
    ComboHitCollection ordChCol;
    ordChCol.reserve(size);
    
    for (int i=0; i<size; ++i) {
      loc = shIndices[i];
      const ComboHit& ch  = (*_hfResult._chcol)[loc];
      if(ch.flag().hasAnyProperty(_hsel) && !ch.flag().hasAnyProperty(_hbkg) ) {
	ordChCol.push_back(ComboHit(ch));
      }
    }
    std::sort(ordChCol.begin(), ordChCol.end(),panelcomp());//zcomp());


    if (_debug>0){
      printf("[RobustHelixFinder::FillHits]-----------------------------------------------------------\n");
      printf("[RobustHelixFinder::FillHits]     i     Face     Panel      X         Y         Z        \n");
      printf("[RobustHelixFinder::FillHits]-----------------------------------------------------------\n");
    }
    
    for (unsigned i=0; i<ordChCol.size(); ++i) {
      // loc = shIndices[i];
      // const ComboHit& ch  = _hfResult._chcol->at(loc);
      ComboHit& ch = ordChCol[i];

      //    if(ch.flag().hasAnyProperty(_hsel) && !ch.flag().hasAnyProperty(_hbkg) ) {
      ComboHit hhit(ch);
      hhit._flag.clear(StrawHitFlag::resolvedphi);
	
      _hfResult._chHitsToProcess.push_back(hhit);

      cx.Station                 = ch.strawId().station();//straw.id().getStation();
      cx.Plane                   = ch.strawId().plane() % 2;//straw.id().getPlane() % 2;
      cx.Face                    = ch.strawId().face();
      cx.Panel                   = ch.strawId().panel();//straw.id().getPanel();

      // get Z-ordered location
      HelixData.orderID(&cx, &co);
     
      int os       = co.Station; 
      int of       = co.Face;
      int op       = co.Panel;

      _hfResult._chHitsWPos.push_back(XYWVec(hhit.pos(),  of, hhit.nStrawHits()));

      int       stationId = os;
      int       faceId    = of + stationId*StrawId::_nfaces*FaceZ_t::kNPlanesPerStation;//RobustHelixFinderData::kNFaces;
      // int       panelId   = op + faceId*RobustHelixDataFinderData::kNPanelsPerFace;
      FaceZ_t* fz        = &HelixData._oTracker[faceId];
      PanelZ_t*pz        = &fz->panelZs[op];

      //	pz->_chHitsToProcess.push_back(hhit);//[fz->fNHits] = hhit;
      //	pz->fNHits  = pz->fNHits + 1;
      if (pz->idChBegin < 0 ){
	pz->idChBegin = _hfResult._chHitsToProcess.size() - 1;
	pz->idChEnd   = _hfResult._chHitsToProcess.size();	
      } else {
	pz->idChEnd   = _hfResult._chHitsToProcess.size();	
      }

      if (fz->idChBegin < 0 ){
	fz->idChBegin = _hfResult._chHitsToProcess.size() - 1;
	fz->idChEnd   = _hfResult._chHitsToProcess.size();	
      } else {
	fz->idChEnd   = _hfResult._chHitsToProcess.size();	
      }
	
      if (_debug>0){
	printf("[RobustHelixFinder::FillHits] %4i %6i %10i %10.3f %10.3f %10.3f\n", nFiltComboHits, faceId, op, ch.pos().x(), ch.pos().y(), ch.pos().z() );
      }
	
      // if (pz->nChHits() > PanelZ_t::kNMaxPanelHits) printf("[RobustHelixDataFinderAlg::fillFaceOrderedHits] number of hits with the panel exceed the limit: NHits =  %i MaxNHits = %i\n", pz->fNHits, PanelZ_t::kNMaxPanelHits);
      ++nFiltComboHits;
      nFiltStrawHits += ch.nStrawHits();
      //      }
    }
    // }
    
    HelixData._nFiltComboHits = nFiltComboHits;  //ComboHit counter
    HelixData._nFiltStrawHits = nFiltStrawHits;  //StrawHit counter
    
    if (_diag) {
      HelixData._diag.nChPPanel = 0;
      HelixData._diag.nChHits   = HelixData._chHitsToProcess.size();

      FaceZ_t*      facez;
      PanelZ_t*     panelz;
    
      int           nhitsFace(0);      

      if (_debug>0){
	printf("[RobustHelixFinder::ReadHits]-----------------------------------------------------------\n");
	printf("[RobustHelixFinder::ReadHits]    i     Face     Panel      X         Y         Z      \n");
	printf("[RobustHelixFinder::ReadHits]-----------------------------------------------------------\n");
	
	for (unsigned i=0; i<HelixData._chHitsToProcess.size(); ++i){
	  ComboHit* ch = &HelixData._chHitsToProcess[i];
	  printf("[RobustHelixFinder::ReadHits] %4i %6i %10i %10.3f %10.3f %10.3f\n", i, ch->strawId().uniqueFace(), ch->strawId().panel(), ch->pos().x(), ch->pos().y(), ch->pos().z() );
	}

	
	printf("[RobustHelixFinder::ReadIndeces]----------------------------------------------------------------------\n");
	printf("[RobustHelixFinder::ReadIndeces]    Face        fBg       fEnd       Panel       pBg        pEnd      \n");
	printf("[RobustHelixFinder::ReadIndeces]----------------------------------------------------------------------\n");

	for (int f=0; f<StrawId::_ntotalfaces; ++f){
	  facez     = &HelixData._oTracker[f];

	  for (int p=0; p<FaceZ_t::kNPanels; ++p){
	    panelz = &facez->panelZs[p];  
	    if (panelz->nChHits() != 0) printf("[RobustHelixFinder::ReadIndeces] %6i %10i %10i %10i %10i %10i\n", f, facez->idChBegin, facez->idChEnd, p, panelz->idChBegin, panelz->idChEnd);
	  }
	}
      }
      for (int f=0; f<StrawId::_ntotalfaces; ++f){
	facez     = &HelixData._oTracker[f];
      
	for (int p=0; p<FaceZ_t::kNPanels; ++p){
	  panelz = &facez->panelZs[p];  
	  nhitsFace = panelz->nChHits();
	  if ( nhitsFace > HelixData._diag.nChPPanel) HelixData._diag.nChPPanel = nhitsFace;
	}//end loop over the panel    
      }//end loop over the faces
    }
  }

  unsigned  RobustHelixFinder::filterCircleHits(RobustHelixFinderData& helixData)
  {
    unsigned changed(0);
    int      nGoodSH(0);
    static XYZVec zaxis(0.0,0.0,1.0); // unit in z direction
    RobustHelix& helix = helixData._hseed._helix;

    // loop over hits
    ComboHit*     hit(0);
    FaceZ_t*      facez;
    
    float         chi2dXY(0), chCounter(1e-10);
    int           nhitsFace(0);

    //for diagnostic purposes
    float         drBestVec   [RobustHelixFinderData::kMaxResidIndex]={-9999.};
    float         rwdotBestVec[RobustHelixFinderData::kMaxResidIndex]={-9999.};

    for (int f=0; f<StrawId::_ntotalfaces; ++f){
      facez     = &helixData._oTracker[f];
      
      float      minChi2(_maxrpull);
      float      drBest(-1.), rwdotBest(-1.); 
      HitInfo_t  indexBestComboHit;
      bool       oldoutBest(false);
      
      nhitsFace = facez->nChHits();
      if (nhitsFace == 0)                        continue;

      for (int ip=0; ip<nhitsFace; ++ip){
	hit = &helixData._chHitsToProcess[facez->idChBegin + ip];

	bool oldout = hit->_flag.hasAnyProperty(_outlier);
	hit->_flag.clear(_outlier);

	const XYZVec& wdir = hit->wdir();
	XYZVec cvec = PerpVector(hit->pos() - helix.center(),Geom::ZDir()); // direction from the circle center to the hit
	XYZVec cdir = cvec.Unit(); // direction from the circle center to the hit
	float rwdot = wdir.Dot(cdir); // compare directions of radius and wire
	if(rwdot > _maxrwdot){
	  hit->_flag.merge(_outlier);
	  //	  if(!oldout) ++changed;
	  continue;
	}
	float dr = sqrtf(cvec.mag2())-helix.radius();
	if ( fabs(dr) > _maxdr ) {
	  hit->_flag.merge(_outlier);
	  //	  if(!oldout) ++changed;
	  continue;
	}
	
	float rwdot2 = rwdot*rwdot;
	// compute radial difference and pull
	float werr = hit->posRes(StrawHitPosition::wire);
	float terr = hit->posRes(StrawHitPosition::trans);
	// the resolution is dominated the resolution along the wire
	//	float rres = std::max(sqrtf(werr*werr*rwdot2 + terr*terr*(1.0-rwdot2)),_minrerr);
	float rres = sqrtf(werr*werr*rwdot2 + terr*terr*(1.0-rwdot2));
	float rpull = fabs(dr/rres)*_rpullScaleF;
	if ( rpull > _maxrpull ) {
	  hit->_flag.merge(_outlier);
	  //	  if(!oldout) ++changed;
	  continue;
	}
	
	//	if (oldout) ++changed;

	if ( rpull < minChi2){
	  minChi2    = rpull;
	  drBest     = dr;
	  rwdotBest  = rwdot; 
	  oldoutBest = oldout;
	  indexBestComboHit.face          = f;
	  indexBestComboHit.panel         = hit->strawId().uniquePanel();
	  indexBestComboHit.panelHitIndex = facez->idChBegin + ip;
	}
	  
	//set all the hits as outlier. Only the best within the face will be cleared
	hit->_flag.merge(_outlier);
	  
      }//end loop over the panels
      
      if (indexBestComboHit.face >=0 ) {
	hit     = &helixData._chHitsToProcess[indexBestComboHit.panelHitIndex];
	
	//remove the outlier flag
	hit->_flag.clear(StrawHitFlag::outlier);
	nGoodSH += hit->nStrawHits();
	chi2dXY += minChi2*minChi2;
	
	if(oldoutBest) ++changed;

	if (chCounter < RobustHelixFinderData::kMaxResidIndex) {
	  drBestVec   [int(chCounter)] = drBest;
	  rwdotBestVec[int(chCounter)] = rwdotBest;
	}
	chCounter += 1.;
      }
      
    }//end loop over the faces
    
    helixData._nXYSh = nGoodSH;
    helixData._nXYCh = int(chCounter);

    helix._chi2dXY   = chi2dXY/chCounter;

    if (_diag) {
      helixData._diag.chi2dXY = chi2dXY/chCounter;
      helixData._diag.nXYCh   = helixData._nXYCh;

      for (int i=0; i<int(chCounter); ++i){
	if (drBestVec   [i]>-999.) helixData._diag.resid[i] = drBestVec   [i];
	if (rwdotBestVec[i]>-999.) helixData._diag.rwdot[i] = rwdotBestVec[i];
      }
      //      helixData._diag. = chi2dXY/chCounter;
    }
    
    return changed;
  }


  void RobustHelixFinder::fitHelix(RobustHelixFinderData& helixData){
    // iteratively fit the helix including filtering
    unsigned niter(0);
    unsigned nitermva(0);
    bool     changed(true), xychanged(true), fzchanged(true);
    unsigned niterxy(0), niterfz(0);

    do {
      niterxy = 0;
      do {
	_hfit.fitCircle(helixData, _targetcon);
	xychanged = filterCircleHits(helixData) > 0;
	++niterxy;
      } while (helixData._hseed._status.hasAllProperties(TrkFitFlag::circleOK) && niterxy < _maxniter && xychanged);
   
      if (_diag) {
	helixData._diag.xyniter  = niterxy;
	helixData._diag.nShFitXY = helixData._nXYSh;
	helixData._diag.nChFitXY = helixData._nXYCh;
      }
      
      if (helixData._nXYSh < _minnsh) {
	helixData._hseed._status.clear(TrkFitFlag::circleOK);
	niter = _maxniter;//exit from this while()
      }

      // then fit phi-Z
      if (helixData._hseed._status.hasAnyProperty(TrkFitFlag::circleOK)) {
	if (niterxy < _maxniter)
	  helixData._hseed._status.merge(TrkFitFlag::circleConverged);
	else
	  helixData._hseed._status.clear(TrkFitFlag::circleConverged);

	// solve for the longitudinal parameters
	niterfz = 0;
	fzchanged = false;
	do {
	  _hfit.fitFZ(helixData);
	  fzchanged = filterHits(helixData);
	  ++niterfz;
	} while (helixData._hseed._status.hasAllProperties(TrkFitFlag::phizOK)  && niterfz < _maxniter && fzchanged);

	if (helixData._nZPhiSh < _minnsh) {
	  helixData._hseed._status.clear(TrkFitFlag::phizOK);
	  niter = _maxniter;//exit from this while()
	}

	if (helixData._hseed._status.hasAnyProperty(TrkFitFlag::phizOK)) {
	  if (niterfz < _maxniter)
	    helixData._hseed._status.merge(TrkFitFlag::phizConverged);
	  else
	    helixData._hseed._status.clear(TrkFitFlag::phizConverged);
	}
      }
      //here is where we should check for the hits within the face to searchfor missing/best ones
      ++niter;
      changed = fzchanged || xychanged;

      // update the stereo hit positions; this checks how much the positions changed
      // do this only in non trigger mode

      if (_updateStereo && _hfit.goodHelix(helixData._hseed.helix()))
	changed |= updateStereo(helixData);
    } while (_hfit.goodHelix(helixData._hseed.helix()) && niter < _maxniter && changed);

    if (_diag) helixData._diag.niter = niter;

    if (_hfit.goodHelix(helixData._hseed.helix())  && 
	helixData._hseed._status.hasAnyProperty(TrkFitFlag::circleOK) && 
	helixData._hseed._status.hasAnyProperty(TrkFitFlag::phizOK) ) {
	
      helixData._hseed._status.merge(TrkFitFlag::helixOK);
      updateT0(helixData);
      if (niter < _maxniter) helixData._hseed._status.merge(TrkFitFlag::helixConverged);

      if (_usemva) {
	bool changed = true;
	while (helixData._hseed._status.hasAllProperties(TrkFitFlag::helixOK)  && nitermva < _maxniter && changed) {
	  fillMVA(helixData);
	  changed = filterHitsMVA(helixData);
	  if (!changed) break;
	  refitHelix(helixData);
	  // update t0 each iteration as that's used in the MVA
	  updateT0(helixData);
	  ++nitermva;
	}
	if (nitermva < _maxniter)
	  helixData._hseed._status.merge(TrkFitFlag::helixConverged);
	else
	  helixData._hseed._status.clear(TrkFitFlag::helixConverged);
      }
    }
    if (_diag > 0){
      _niter->Fill(niter);
      _niterfz->Fill(niterfz);
      _niterxy->Fill(niterxy);
      _nitermva->Fill(nitermva);
      if (!_usemva) fillMVA(helixData);
    }
  }



  //------------------------------------------------------------------------------------------


  void RobustHelixFinder::fitChi2Helix(RobustHelixFinderData& helixData){
    // iteratively fit the helix including filtering

    //before starting, try to resolve the z-phi part of the helix 
    _chi2hfit.initFitChi2FZ(helixData);

    //use chi2 line fit to make some preliminary cleanup
    _chi2hfit.fitChi2FZ(helixData,0);
    _chi2hfit.fitChi2FZ(helixData);

    unsigned xyniter(0);
    int      xychanged = filterChi2XYHits(helixData);
    while (helixData._hseed._status.hasAnyProperty(TrkFitFlag::circleOK) && xyniter < _maxniter && (xychanged!=0)) {
      xychanged = filterChi2XYHits(helixData);
      ++xyniter;
    } 
   
    if (_diag) {
      helixData._diag.xyniter = xyniter;
      helixData._diag.nShFitXY   = helixData._nXYSh;
      helixData._diag.nChFitXY   = helixData._nXYCh;
    }

    // then fit phi-Z
    if (helixData._hseed._status.hasAnyProperty(TrkFitFlag::circleOK)) {
      if (xyniter < _maxniter)
	helixData._hseed._status.merge(TrkFitFlag::circleConverged);
      else
	helixData._hseed._status.clear(TrkFitFlag::circleConverged);

      // solve for the longitudinal parameters
      unsigned fzniter(0);
      _chi2hfit.fitChi2FZ(helixData);
      int fzchanged = filterChi2ZPhiHits(helixData);
      while (helixData._hseed._status.hasAnyProperty(TrkFitFlag::phizOK)  && fzniter < _maxniter && (fzchanged!=0)) {
	fzchanged = filterChi2ZPhiHits(helixData);
	++fzniter;
      }

      if (_diag) helixData._diag.fzniter = fzniter;

      if (helixData._hseed._status.hasAnyProperty(TrkFitFlag::phizOK)) {
	if (fzniter < _maxniter){
	  helixData._hseed._status.merge(TrkFitFlag::phizConverged);
	  
	  //now update all the helix parameters
	  updateChi2HelixInfo(helixData);

	  if (_hfit.goodHelix(helixData._hseed.helix()) && _chi2hfit.goodHelixChi2(helixData)) {
	    helixData._hseed._status.merge(TrkFitFlag::helixOK);

	    //now search for missing hits
	    findMissingHits(helixData);

	    updateT0(helixData);
      
	    helixData._hseed._status.merge(TrkFitFlag::helixConverged);

	    _chi2hfit.defineHelixParams(helixData);
	  }
	}
	else
	  helixData._hseed._status.clear(TrkFitFlag::phizConverged);
      }
      
    }

  }

  void RobustHelixFinder::refitHelix(RobustHelixFinderData& helixData) {
    // reset the fit status flags, in case this is called iteratively
    helixData._hseed._status.clear(TrkFitFlag::helixOK);      
    _hfit.fitCircle(helixData, _targetcon);
    if (helixData._hseed._status.hasAnyProperty(TrkFitFlag::circleOK)) {
      _hfit.fitFZ(helixData);
      if (_hfit.goodHelix(helixData._hseed._helix)) helixData._hseed._status.merge(TrkFitFlag::helixOK);
    }
  }

  unsigned RobustHelixFinder::hitCount(RobustHelixFinderData& helixData) {
    unsigned nHits(0);
    
    ComboHit*      hit(0);
    
    for (unsigned f=0; f<helixData._chHitsToProcess.size(); ++f){
      hit = &helixData._chHitsToProcess[f];
      if (hit->_flag.hasAnyProperty(_outlier))   continue;
      ++nHits;
    }//end faces loop

    return nHits;
  }


  bool RobustHelixFinder::updateStereo(RobustHelixFinderData& helixData) {
    static StrawHitFlag stereo(StrawHitFlag::stereo);
    bool retval(false);
    // loop over the stereo hits in the helix and update their positions given the local helix direction
    for(auto& ch : helixData._hseed._hhits){
      if(ch.flag().hasAllProperties(stereo) && ch.nCombo() >=2) {
	// local helix direction at the average z of this hit
	XYZVec hdir;
	helixData._hseed.helix().direction(ch.pos().z(),hdir);
	// needs re-implementing with ComboHits FIXME!	  
	//	XYZVec pos1, pos2;
	//	sthit.position(shcol,tracker,pos1,pos2,hdir);
      }
    }
    return retval;
  }


  void RobustHelixFinder::fillPluginDiag(RobustHelixFinderData& helixData, int helCounter) {
    //--------------------------------------------------------------------------------    
    // fill diagnostic information
    //--------------------------------------------------------------------------------
    float          mm2MeV = 3./10.;  // approximately , at B=1T
	  
    int loc = _data.nseeds[helCounter];
    if (loc < _data.maxSeeds()) {
      _data.nChPPanel   [helCounter][loc] = helixData._diag.nChPPanel;
      _data.nChHits     [helCounter][loc] = helixData._diag.nChHits;

      int nhits          = helixData._hseed._hhits.size();
      _data.ntclhits    [helCounter][loc] = helixData._timeCluster->hits().size();
      _data.nhits       [helCounter][loc] = nhits;

      _data.ntriplet0   [helCounter][loc] = helixData._diag.ntriple_0;
      _data.ntriplet1   [helCounter][loc] = helixData._diag.ntriple_1;
      _data.ntriplet2   [helCounter][loc] = helixData._diag.ntriple_2;

      _data.xyniter     [helCounter][loc] = helixData._diag.xyniter;
      _data.fzniter     [helCounter][loc] = helixData._diag.fzniter;
      _data.niter       [helCounter][loc] = helixData._diag.niter;
      _data.nrescuedhits[helCounter][loc] = helixData._diag.nrescuedhits;
   
      _data.nShFitCircle[helCounter][loc] = helixData._diag.nShFitCircle;
      _data.nChFitCircle[helCounter][loc] = helixData._diag.nChFitCircle;
      _data.nShFitXY    [helCounter][loc] = helixData._diag.nShFitXY;
      _data.nChFitXY    [helCounter][loc] = helixData._diag.nChFitXY;


      _data.nfz0counter [helCounter][loc] = helixData._diag.nfz0counter;

      _data.nshsxy_0    [helCounter][loc] = helixData._diag.nshsxy_0;
      _data.rsxy_0      [helCounter][loc] = helixData._diag.rsxy_0;
      _data.chi2dsxy_0  [helCounter][loc] = helixData._diag.chi2dsxy_0;
   	                
      _data.nshsxy_1    [helCounter][loc] = helixData._diag.nshsxy_1;
      _data.rsxy_1      [helCounter][loc] = helixData._diag.rsxy_1;
      _data.chi2dsxy_1  [helCounter][loc] = helixData._diag.chi2dsxy_1;
   	                
      _data.nshszphi_0  [helCounter][loc] = helixData._diag.nshszphi_0;
      _data.lambdaszphi_0    [helCounter][loc] = helixData._diag.lambdaszphi_0;
      _data.chi2dszphi_0[helCounter][loc] = helixData._diag.chi2dszphi_0;
   	                
      _data.nshszphi_1  [helCounter][loc] = helixData._diag.nshszphi_1;
      _data.lambdaszphi_1    [helCounter][loc] = helixData._diag.lambdaszphi_1;
      _data.chi2dszphi_1[helCounter][loc] = helixData._diag.chi2dszphi_1;

   
      _data.nXYSh       [helCounter][loc] = helixData._nXYSh;
      _data.nZPhiSh     [helCounter][loc] = helixData._nZPhiSh;

      _data.rinit       [helCounter][loc] = helixData._diag.radius_0;
      _data.lambda0     [helCounter][loc] = helixData._diag.lambda_0;
      _data.lambda1     [helCounter][loc] = helixData._diag.lambda_1;
      _data.radius      [helCounter][loc] = helixData._hseed.helix().radius();
      _data.pT          [helCounter][loc] = mm2MeV*_data.radius[helCounter][loc];
      _data.p           [helCounter][loc] = _data.pT[helCounter][loc]/std::cos( std::atan(helixData._hseed.helix().lambda()/_data.radius[helCounter][loc]));
	
      _data.chi2XY      [helCounter][loc] = helixData._diag.chi2dXY;
      _data.chi2ZPhi    [helCounter][loc] = helixData._diag.chi2dZPhi;
	    
      _data.nseeds[helCounter]++;
	    
      _data.dr          [helCounter][loc] = helixData._diag.radius_2 - helixData._diag.radius_1;
      _data.chi2d_helix [helCounter][loc] = helixData._diag.chi2d_helix;
      
      _data.nXYCh       [helCounter][loc] = helixData._diag.nXYCh;

      _data.nLoops      [helCounter][loc] = helixData._diag.nLoops           ;

      _data.nHitsLoopFailed[helCounter][loc] = helixData._diag.nHitsLoopFailed;

      _data.meanHitRadialDist [helCounter][loc] = helixData._diag.meanHitRadialDist;

      for (int i=0; i<helixData._diag.nXYCh; ++i) {
	if (helixData._diag.rwdot[i]>-999.) _data.hitRWDot[helCounter][loc][i] = helixData._diag.rwdot[i];
      	
	if (helixData._diag.resid[i]>-999.) _data.hitDr   [helCounter][loc][i] = helixData._diag.resid[i];
	else break;
      }
    }   else {
      printf(" N(seeds) > %i, IGNORE SEED\n",_data.maxSeeds());
    }
  }

  void     RobustHelixFinder::updateChi2HelixInfo(RobustHelixFinderData& helixData){
    RobustHelix&  helix        = helixData._hseed._helix;
    XYVec         center       = XYVec(helix.center().x(), helix.center().y());
    float         radius       = helixData._sxy.radius();
    center.SetX(helixData._sxy.x0());
    center.SetY(helixData._sxy.y0());

    //update the LSqsums
    ComboHit*      hit(0);
    
    helixData._sxy.clear();
    if (_targetcon) helixData._sxy.addPoint(0.,0.,1./900.);

    helixData._szphi.clear();
    helixData._nXYSh   = 0;
    helixData._nZPhiSh = 0;

    for (unsigned f=0; f<helixData._chHitsToProcess.size(); ++f){
      hit    = &helixData._chHitsToProcess[f];
	 
      if (hit->_flag.hasAnyProperty(_outlier))   continue;
	
      hit->_xyWeight   = _chi2hfit.evalWeightXY(*hit, center);
      hit->_zphiWeight = _chi2hfit.evalWeightZPhi(*hit, center, radius);
	  
      helixData._sxy.addPoint(hit->pos().x(), hit->pos().y(), hit->_xyWeight);
      helixData._szphi.addPoint(hit->pos().z(), hit->helixPhi(), hit->_zphiWeight);

      helixData._nXYSh   += hit->nStrawHits();
      helixData._nZPhiSh += hit->nStrawHits();
    }
    
    center.SetX(helixData._sxy.x0());
    center.SetY(helixData._sxy.y0());

    //update the circle part
    helix._rcent  = sqrtf(center.Mag2());
    helix._fcent  = polyAtan2(center.y(), center.x());
    helix._radius = helixData._sxy.radius();
 
    //update the Z-Phi part
    helix._lambda = 1./(helixData._szphi.dfdz());
    helix._fz0    = helixData._szphi.phi0();
  }
  
  void     RobustHelixFinder::updateHelixXYInfo(RobustHelixFinderData& helixData){
    RobustHelix&  helix        = helixData._hseed._helix;
    XYVec         center       = XYVec(helix.center().x(), helix.center().y());

    center.SetX(helixData._sxy.x0());
    center.SetY(helixData._sxy.y0());

    helix._rcent  = sqrtf(center.Mag2());
    helix._fcent  = polyAtan2(center.y(), center.x());
    helix._radius = helixData._sxy.radius();
  }

  void     RobustHelixFinder::updateHelixZPhiInfo(RobustHelixFinderData& helixData){
    RobustHelix&  helix        = helixData._hseed._helix;

    helix._lambda = 1./(helixData._szphi.dfdz());
    helix._fz0    = helixData._szphi.phi0();
  }

  //----------------------------------------------------------------------------------------------------
  void     RobustHelixFinder::searchWorstHitXY(RobustHelixFinderData& helixData, HitInfo_t& hitInfo){
    hitInfo.face          = -1;
    hitInfo.panel         = -1;
    hitInfo.panelHitIndex = -1;
    
    float          dr, wt, hitChi2;
    XYVec          centerXY     = XYVec(helixData._hseed._helix.center().x(), helixData._hseed._helix.center().y());
    XYZVec         center       = helixData._hseed._helix.center();
    float          helix_radius = helixData._hseed._helix.radius();
    float          hitChi2Worst = _maxrpull*_maxrpull;

    ComboHit*      hit(0);
    FaceZ_t*       facez(0);
    int            nhitsFace(0);
    
    for (int f=0; f<StrawId::_ntotalfaces; ++f){
      facez     = &helixData._oTracker[f];
      nhitsFace = facez->nChHits();
      
      if (nhitsFace == 0)                          continue;

      for (int ip=0; ip<nhitsFace; ++ip){
	hit = &helixData._chHitsToProcess[facez->idChBegin + ip];
	
	if (hit->_flag.hasAnyProperty(_outlier))   continue;
	
	XYZVec  cvec  = PerpVector(hit->pos() - center,Geom::ZDir());
	dr    = sqrtf(cvec.mag2()) - helix_radius;
	wt    = _chi2hfit.evalWeightXY(*hit, centerXY);

	hitChi2 = dr*dr*wt;

	if (hitChi2 > hitChi2Worst) {
	  hitChi2Worst          = hitChi2;
	  hitInfo.face          = f;
	  hitInfo.panel         = hit->strawId().uniquePanel();
	  hitInfo.panelHitIndex = facez->idChBegin + ip;
	  hit->_xyWeight        = wt;
	}
      }//end loop of the hits within the face
    }//end faces loop
  }
  
  //----------------------------------------------------------------------------------------------------
  void     RobustHelixFinder::searchWorstHitZPhi(RobustHelixFinderData& helixData, HitInfo_t& hitInfo){
    hitInfo.face          = -1;
    hitInfo.panel         = -1;
    hitInfo.panelHitIndex = -1;
    
    // XYZVec         center       = helixData._hseed._helix.center();
    XYVec          centerXY     = XYVec(helixData._hseed._helix.center().x(), helixData._hseed._helix.center().y());
    float          helix_radius = helixData._hseed._helix.radius();

    ComboHit*      hit(0);
    FaceZ_t*       facez(0);
    int            nhitsFace(0);
    
    float          chi2min(1e10), chi2;
    ::LsqSums4     szphi;

    for (int f=0; f<StrawId::_ntotalfaces; ++f){
      facez     = &helixData._oTracker[f];
      nhitsFace = facez->nChHits();

      if (nhitsFace == 0)                          continue;

      for (int ip=0; ip<nhitsFace; ++ip){
	hit = &helixData._chHitsToProcess[facez->idChBegin + ip];
	
	if (hit->_flag.hasAnyProperty(_outlier))   continue;
	
	szphi.init(helixData._szphi);

	// XYZVec  cvec  = PerpVector(hit->pos() - center,Geom::ZDir());
	float   phi   = hit->helixPhi();
	float   wt    = _chi2hfit.evalWeightZPhi(*hit, centerXY, helix_radius);

	szphi.removePoint(hit->pos().z(), phi, wt);
	chi2 = szphi.chi2DofLine();

	if (chi2 < chi2min) {
	  chi2min               = chi2;
	  hitInfo.face          = f;
	  hitInfo.panel         = hit->strawId().uniquePanel();
	  hitInfo.panelHitIndex = facez->idChBegin + ip;

	  hit->_zphiWeight      = wt;
	  //	    hitInfo.weightZPhi    = wt;
	}
      }//end loop pver the hits within a face
    }//end faces loop
  }
}
using mu2e::RobustHelixFinder;
DEFINE_ART_MODULE(RobustHelixFinder);
