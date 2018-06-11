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
#include "GeometryService/inc/GeomHandle.hh"
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

#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Mu2eUtilities/inc/ModuleHistToolBase.hh"
#include "Mu2eUtilities/inc/polyAtan2.hh"
#include "art/Utilities/make_tool.h"

#include "TrkPatRec/inc/RobustHelixFinder_types.hh"
#include "TrkReco/inc/RobustHelixFinderData.hh"

#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"

#include <boost/accumulators/accumulators.hpp>
#include "boost_fix/accumulators/statistics/stats.hpp"
#include "boost_fix/accumulators/statistics.hpp"
#include <boost/accumulators/statistics/median.hpp>

#include "TH1F.h"
#include "Math/VectorUtil.h"

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
    int                                 _diag,_debug;
    int                                 _printfreq;
    bool				_prefilter; // prefilter hits based on sector
    bool				_updatestereo; // update the stereo hit positions each iteration
    int 				_minnsh; // minimum # of strawHits to work with
    unsigned				_minnhit; // minimum # of hits to work with
    float                               _maxchi2dxy;
    float                               _maxchi2dzphi;
    float                               _maxphihitchi2;
    float				_maxdr; // maximum hit-helix radius difference
    float				_maxrpull; // maximum hit-helix radius difference pull
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
    float				_minmva; // outlier cut on MVA

    art::InputTag			_ccTag;
    art::InputTag			_chTag;
    art::InputTag			_chfTag;
    art::InputTag			_tcTag;

    StrawHitFlag  _hsel, _hbkg;

    MVATools _stmva, _nsmva;
    HelixHitMVA _vmva; // input variables to TMVA for filtering hits

    TH1F* _niter, *_nitermva;

    RobustHelixFit   _hfit;
    std::vector<Helicity> _hels; // helicity values to fit 
    TrkTimeCalculator _ttcalc;
    float             _t0shift;   
    StrawHitFlag      _outlier;
    bool              _updateStereo;
    
    std::unique_ptr<ModuleHistToolBase>   _hmanager;
    RobustHelixFinderTypes::Data_t        _data;
    RobustHelixFinderData                 _hfResult;


    void     findHelices(ComboHitCollection& chcol, const TimeClusterCollection& tccol);    
    void     prefilterHits(RobustHelixFinderData& helixData, int& nFilteredStrawHits); 
    unsigned filterCircleHits(RobustHelixFinderData& helixData); 
    int      filterZPhiHits(RobustHelixFinderData& helixData);
    int      filterXYHits(RobustHelixFinderData& helixData); 
    bool     filterHits(RobustHelixFinderData& helixData);
    void     fillMVA(RobustHelixFinderData& helixData); 
    bool     filterHitsMVA(RobustHelixFinderData& helixData);
    void     updateT0(RobustHelixFinderData& helixData);
    bool     updateStereo(RobustHelixFinderData& helixData);
    unsigned hitCount(RobustHelixFinderData& helixData);
    void     fillFaceOrderedHits(RobustHelixFinderData& helixData);
    void     fillGoodHits(RobustHelixFinderData& helixData);
    void     fitHelix(RobustHelixFinderData& helixData);
    void     fitHelix_2(RobustHelixFinderData& helixData);
    void     refitHelix(RobustHelixFinderData& helixData);
    void     findMissingHits(RobustHelixFinderData& helixData);
    void     fillPluginDiag(RobustHelixFinderData& helixData, int helCounter);
    void     updateHelixInfo    (RobustHelixFinderData& helixData);
    void     updateHelixXYInfo  (RobustHelixFinderData& helixData);
    void     updateHelixZPhiInfo(RobustHelixFinderData& helixData);
    void     searchWorstHitXY   (RobustHelixFinderData& helixData, RobustHelixFinderData::HitInfo_t& hitInfo);
    void     searchWorstHitZPhi (RobustHelixFinderData& helixData, RobustHelixFinderData::HitInfo_t& hitInfo);
  };

  RobustHelixFinder::RobustHelixFinder(fhicl::ParameterSet const& pset) :
    _diag        (pset.get<int>("diagLevel",0)),
    _debug       (pset.get<int>("debugLevel",0)),
    _printfreq   (pset.get<int>("printFrequency",101)),
    _prefilter   (pset.get<bool>("PrefilterHits",true)),
    _updatestereo(pset.get<bool>("UpdateStereoHits",false)),
    _minnsh      (pset.get<int>("minNStrawHits",10)),
    _minnhit	 (pset.get<unsigned>("minNHit",5)),
    _maxchi2dxy  (pset.get<float>("MaxChi2dXY", 5.0)),
    _maxchi2dzphi(pset.get<float>("MaxChi2dZPhi", 5.0)),
    _maxphihitchi2(pset.get<float>("MaxHitPhiChi2", 25.0)),
    _maxdr	 (pset.get<float>("MaxRadiusDiff",100.0)), // mm
    _maxrpull	 (pset.get<float>("MaxRPull",5.0)), // unitless
    _maxphisep	 (pset.get<float>("MaxPhiHitSeparation",1.0)),
    _saveflag    (pset.get<vector<string> >("SaveHelixFlag",vector<string>{"HelixOK"})),
    _maxniter    (pset.get<unsigned>("MaxIterations",10)), // iterations over outlier removal
    _cradres	 (pset.get<float>("CenterRadialResolution",20.0)),
    _cperpres	 (pset.get<float>("CenterPerpResolution",12.0)),
    _maxdwire    (pset.get<float>("MaxWireDistance",200.0)), // max distance along wire
    _maxdtrans   (pset.get<float>("MaxTransDistance",80.0)), // max distance perp to wire (and z)
    _maxchisq    (pset.get<float>("MaxChisquared",100.0)), // max chisquared
    _maxrwdot	(pset.get<float>("MaxRWDot",1.0)),
    _minrerr     (pset.get<float>("MinRadiusErr",20.0)), // mm
    _usemva      (pset.get<bool>("UseHitMVA",false)),
    _minmva      (pset.get<float> ("MinMVA",0.1)), // min MVA output to define an outlier
    _ccTag	 (pset.get<art::InputTag>("CaloClusterCollection","CaloClusterFast")),
    _chTag	 (pset.get<art::InputTag>("ComboHitCollection")),
    _chfTag	 (pset.get<art::InputTag>("ComboHitFlagCollection")),
    _tcTag	 (pset.get<art::InputTag>("TimeClusterCollection")),
    _hsel        (pset.get<std::vector<std::string> >("HitSelectionBits",std::vector<string>{"TimeDivision"})),
    _hbkg        (pset.get<std::vector<std::string> >("HitBackgroundBits",std::vector<std::string>{"Background"})),
    _stmva       (pset.get<fhicl::ParameterSet>("HelixStereoHitMVA",fhicl::ParameterSet())),
    _nsmva       (pset.get<fhicl::ParameterSet>("HelixNonStereoHitMVA",fhicl::ParameterSet())),
    _hfit        (pset.get<fhicl::ParameterSet>("RobustHelixFit",fhicl::ParameterSet())),
    _ttcalc      (pset.get<fhicl::ParameterSet>("T0Calculator",fhicl::ParameterSet())),
    _t0shift     (pset.get<float>("T0Shift",4.0)),
    _outlier     (StrawHitFlag::outlier),
    _updateStereo    (pset.get<bool>("UpdateStereo",true))
  {
    std::vector<int> helvals = pset.get<std::vector<int> >("Helicities",vector<int>{Helicity::neghel,Helicity::poshel});
    for(auto hv : helvals) {
      Helicity hel(hv);
      _hels.push_back(hel);
      produces<HelixSeedCollection>(Helicity::name(hel));

      if (_diag != 0) _hmanager = art::make_tool<ModuleHistToolBase>(pset.get<fhicl::ParameterSet>("diagPlugin"));
      else            _hmanager = std::make_unique<ModuleHistToolBase>();
    }

    _data.result    = &_hfit;

  }
  
  RobustHelixFinder::~RobustHelixFinder(){}

  //-----------------------------------------------------------------------------
  void RobustHelixFinder::beginRun(art::Run& ) {
    mu2e::GeomHandle<mu2e::TTracker> th;
    const TTracker* tracker = th.get();

    // mu2e::GeomHandle<mu2e::Calorimeter> ch;
    // _calorimeter = ch.get();
    // 					// calibrations

    // mu2e::ConditionsHandle<TrackerCalibrations> tcal("ignored");
    // _trackerCalib = tcal.operator ->();

    // _hfinder.setTracker    (_tracker);
    // _hfinder.setCalorimeter(_calorimeter);

    mu2e::RobustHelixFinderData::ChannelID cx, co;

    for (int ist=0; ist<tracker->nStations(); ist++) {
      const Station* st = &tracker->getStation(ist);
      
      for (int ipl=0; ipl<st->nPlanes(); ipl++) {
	const Plane* pln = &st->getPlane(ipl);
	for (int ipn=0; ipn<pln->nPanels(); ipn++) {
	  const Panel* panel = &pln->getPanel(ipn);
	  int face;
	  if (panel->id().getPanel() % 2 == 0) face = 0;
	  else                                 face = 1;
	  cx.Station = ist;
	  cx.Plane   = ipl;
	  cx.Face    = face;
	  cx.Panel   = ipn;
	  //	    cx.Layer   = il;
	  _hfResult.orderID (&cx, &co);
	  int os = co.Station; 
	  int of = co.Face;
	  // int op = co.Panel;

	  int       stationId = os;
	  int       faceId    = of + stationId*RobustHelixFinderData::kNFaces;
	  // int       panelId   = op + faceId*CalHelixFinderData::kNPanelsPerFace;
	  RobustHelixFinderData::FaceZ_t*  fz        = &_hfResult._oTracker[faceId];

	  //-----------------------------------------------------------------------------
	  // face caches the z coordinate
	  //-----------------------------------------------------------------------------
	  fz->z      = (panel->getStraw(0).getMidPoint().z()+panel->getStraw(1).getMidPoint().z())/2.;
	  fz->fNHits = 0;
	}	
      }
    }

    if (_debug > 10){
      printf("//-------------------------//\n");
      printf("//     Face      Z        //\n");
      printf("//-------------------------//\n");

      RobustHelixFinderData::FaceZ_t* facez(0);
    
      for (int p=0; p<RobustHelixFinderData::kNTotalPanels; ++p){
	facez = &_hfResult._oTracker[p];
	double z = facez->z;
	printf("//  %5i     %10.3f //\n", p, z);
      }
      printf("//----------------------------------//\n");

    }
    
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
      _nitermva = tfs->make<TH1F>( "nitermva" , "Number of MVA Fit Iteraions",201,-0.5,200.5);
      _hmanager->bookHistograms(tfs);
    }
  }

  void RobustHelixFinder::produce(art::Event& event ) {
    // find input
    auto tcH = event.getValidHandle<TimeClusterCollection>(_tcTag);
    const TimeClusterCollection& tccol(*tcH);

    art::Handle<ComboHitCollection> chH;
    if(!event.getByLabel(_chTag, chH))
      throw cet::exception("RECO")<<"RobustHelixFinder: No ComboHit collection found for tag" <<  _chTag << endl;
    const ComboHitCollection& chcol(*chH);

    auto chfH = event.getValidHandle<StrawHitFlagCollection>(_chfTag);
    
    const StrawHitFlagCollection*       _chfcol = chfH.product();
    // if(!event.getByLabel(_chfTag, chfH)){
    //   _chfcol= chfH.product();
    // }else {
    //   throw cet::exception("RECO")<<"RobustHelixFinder: No StrawHitFlag collection found for tag" <<  _chfTag << endl;
    // }

    // create output: seperate by helicity
    std::map<Helicity,unique_ptr<HelixSeedCollection>> helcols;
    int counter(0);
    for( auto const& hel : _hels) {
      helcols[hel] = unique_ptr<HelixSeedCollection>(new HelixSeedCollection());
      //      _data.helices[counter] = helcols[hel].get();
      _data.nseeds [counter] = 0;
      ++counter;
    }
    
    _data.event       = &event;
    _data.result      = &_hfit;
    _data.nTimePeaks  = tccol.size();

    _hfResult._chcol  = &chcol;
    _hfResult._chfcol = _chfcol;
      
    // create initial helicies from time clusters: to begin, don't specificy helicity
    for (size_t index=0;index< tccol.size();++index) {
      const auto& tclust = tccol[index];
      HelixSeed hseed;

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
      _hfit.fitCircle(_hfResult);

      //check the number of points associated with the result of the circle fit
      // if (_hfResult._nXYSh < _minnsh)                           continue;
      
      if (_hfResult._hseed._status.hasAnyProperty(TrkFitFlag::circleOK)) {
	// loop over helicities. 
	int     helCounter(0);
	for(auto const& hel : _hels ) {
	  // tentatively put a copy with the specified helicity in the appropriate output vector
	  RobustHelixFinderData tmpResult(_hfResult);
	  tmpResult._hseed._helix._helicity = hel;

	  //fit the helix: refine the XY-circle fit + performs the ZPhi fit
	  // it also performs a clean-up of the hits with large residuals
	  fitHelix(tmpResult);
	  // fitHelix_2(tmpResult);

	  if (tmpResult._hseed.status().hasAnyProperty(_saveflag)){
	    //fill the hits in the HelixSeedCollection
	    fillGoodHits(tmpResult);

	    HelixSeedCollection* hcol = helcols[hel].get();
	    hcol->push_back(tmpResult._hseed);
	    if (_diag > 0) {
	      fillPluginDiag(tmpResult, helCounter);
	    }
	  }
	  ++helCounter;
	}
      }	
      
    }
    // put final collections into event 
    if (_diag > 0) _hmanager->fillHistograms(&_data);

    for(auto const& hel : _hels ) {
      event.put(std::move(helcols[hel]),Helicity::name(hel));
    }
  }

  void RobustHelixFinder::fillGoodHits(RobustHelixFinderData& helixData){
    
    RobustHelixFinderData::FaceZ_t*      facez(0);
    ComboHit*                            hit(0);

    RobustHelix& rhel          = helixData._hseed._helix;
    double       z_start(0);
    double       dfdz          = 1./rhel.lambda();
    bool         isFirst(true);

    for (int f=0; f<RobustHelixFinderData::kNTotalFaces; ++f){
      facez = &helixData._oTracker[f];
      int  nhits = facez->fNHits;
      for (int i=0; i<nhits; ++i){   
	hit = &facez->fHitData.at(i);
	// if (hit->_flag.hasAnyProperty(_outlier))     continue;//FIX ME! 
	
	double   hit_z  = hit->pos().z();
	if ( isFirst ){ 
	  z_start = hit_z;
	  isFirst = false;
	}
      
	double   dx     = (hit->pos().x() - rhel.center().y());
	double   dy     = (hit->pos().y() - rhel.center().y());
	double   shphi  = polyAtan2(dy,dx);//XYZVec(hit->pos() - rhel.center()).phi();
	int      nLoops = (hit_z - z_start)/(2.*M_PI/dfdz);
	shphi = shphi + double(nLoops)*2.*M_PI;

	ComboHit                hhit(*hit);
	hhit._hphi = shphi;
	hhit._flag.merge(StrawHitFlag::resolvedphi);
					
	helixData._hseed._hhits.push_back(hhit);
	
      }//end loop over the hits within a face
    }//end loop over the faces
  }


  void RobustHelixFinder::fillMVA(RobustHelixFinderData& helixData)
  {
    RobustHelix& helix = helixData._hseed._helix;

    static XYZVec zaxis(0.0,0.0,1.0); // unit in z direction
    ComboHit*      hhit(0);
    RobustHelixFinderData::FaceZ_t*       facez(0);
    int            nhitsFace(0);
    
    for (int f=0; f<RobustHelixFinderData::kNTotalFaces; ++f){
      facez     = &helixData._oTracker[f];
      nhitsFace = facez->fNHits;

      if (nhitsFace == 0)                          continue;

      for (int ip=0; ip<nhitsFace; ++ip){
	hhit = &facez->fHitData.at(ip);
	
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
	}//end loop over the hits
      }
    }//end loop over the faces
  }

  bool RobustHelixFinder::filterHitsMVA(RobustHelixFinderData& helixData)
  {  
    bool           changed(false);
    ComboHit*      hhit(0);
    RobustHelixFinderData::FaceZ_t*       facez(0);
    int            nhitsFace(0);
    
    for (int f=0; f<RobustHelixFinderData::kNTotalFaces; ++f){
      facez     = &helixData._oTracker[f];
      nhitsFace = facez->fNHits;

      if (nhitsFace == 0)                          continue;

      for (int ip=0; ip<nhitsFace; ++ip){
	hhit = &facez->fHitData.at(ip);
	
	bool oldout = hhit->_flag.hasAnyProperty(_outlier);

	if (hhit->_qual < _minmva ) hhit->_flag.merge(_outlier);
	else                        hhit->_flag.clear(_outlier);

	changed |= oldout != hhit->_flag.hasAnyProperty(_outlier);
      }    
    }//end loop over the faces

    return changed;
  }

  int  RobustHelixFinder::filterXYHits(RobustHelixFinderData& helixData)
  {
    //reset the value of the XY fit result
    helixData._hseed._status.clear(TrkFitFlag::circleOK);

    int           changed(0);
    int           oldNHitsSh = helixData._nXYSh;
    
    RobustHelix&  helix      = helixData._hseed._helix;

    ComboHit*     hit(0);
    RobustHelixFinderData::FaceZ_t*      facez;

    //perform a reduced chi2 fit
    _hfit.refineFitXY(helixData);


    // helixData._hseed._status.clear(TrkFitFlag::circleOK);

    if (helixData._nXYSh >= _minnsh) {//update the helix info
      //need to update the weights in the LSqsum
      _hfit.refineFitXY(helixData);
   
      //      updateHelixXYInfo(helixData);//should be unnecessary!FIXME!
          
      //search and remove the worst hit(s) if necessary
      RobustHelixFinderData::HitInfo_t     worstHit(0,0);
      float      chi2d = helixData._sxy.chi2DofCircle();
      
      while( (chi2d > _maxchi2dxy) && (helixData._nXYSh >= _minnsh) && (worstHit.face >=0)){
	//reset the content of the worstHit
	worstHit.face         = -1;
	worstHit.faceHitIndex = -1;
	worstHit.weightXY     =  0;
	searchWorstHitXY(helixData, worstHit);
	
	if (worstHit.face >=0){//check if a bad was found or not
	  facez = &helixData._oTracker[worstHit.face];
	  hit   = &facez->fHitData.at(worstHit.faceHitIndex);
	  
	  hit->_flag.merge(_outlier);

	  helixData._sxy.removePoint(hit->pos().x(), hit->pos().y(), worstHit.weightXY);
	  helixData._nXYSh -= hit->nStrawHits();
	  _hfit.refineFitXY(helixData);//should be unnecessary!FIXME!
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
  int RobustHelixFinder::filterZPhiHits(RobustHelixFinderData& helixData)
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

    RobustHelixFinderData::FaceZ_t*      facez;

    //    helixData._hseed._status.clear(TrkFitFlag::circleOK);
    _hfit.refineFitZPhi(helixData);

    if (helixData._nZPhiSh >= _minnsh) {//update the helix info
      //need to update the weights in the LSqsum
      _hfit.refineFitZPhi(helixData);

      //      updateHelixZPhiInfo(helixData);//should be unnecessary!FIXME!

      //search and remove the worst hit(s) if necessary
      RobustHelixFinderData::HitInfo_t     worstHit(0,0);
      float      chi2d = helixData._szphi.chi2DofLine();
      
      while( (chi2d > _maxchi2dzphi) && (helixData._nZPhiSh >= _minnsh) && (worstHit.face >=0)){
	//reset the content of the worstHit
	worstHit.face         = -1;
	worstHit.faceHitIndex = -1;
	
	searchWorstHitZPhi(helixData, worstHit);
	
	if (worstHit.face >=0){//check if a bad was found or not
	  facez = &helixData._oTracker[worstHit.face];
	  hit   = &facez->fHitData.at(worstHit.faceHitIndex);
	  
	  hit->_flag.merge(_outlier);

	  helixData._szphi.removePoint(facez->z, hit->helixPhi(), worstHit.weightZPhi);
	  helixData._nZPhiSh -= hit->nStrawHits();
	  _hfit.refineFitZPhi(helixData);
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
    RobustHelixFinderData::FaceZ_t*      facez;
    int           nhitsFace(0);

    for (int f=0; f<RobustHelixFinderData::kNTotalFaces; ++f){
      facez     = &helixData._oTracker[f];
      nhitsFace = facez->fNHits;
      if (nhitsFace == 0)                        continue;

      for (int ip=0; ip<nhitsFace; ++ip){
	hit = &facez->fHitData.at(ip);
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
	    hit->_flag.merge(_outlier);
	    changed = true;
	  }else 
	  {
	    hit->_flag.clear(_outlier);
	    nGoodSH += hit->nStrawHits();
	  }
      }
    }
    
    helixData._nZPhiSh = nGoodSH;
    
    return changed;
  }
    
  void     RobustHelixFinder::findMissingHits(RobustHelixFinderData& helixData){
    RobustHelixFinderData::FaceZ_t*      facez;
    ComboHit*  hit(0);
    int        nhitsPerFace(0), n_added_points(0);
    RobustHelixFinderData::HitInfo_t  bestHit;
    
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
    bestHit.faceHitIndex  = -1;
    bestHit.weightXY      = 1.;
    bestHit.weightZPhi    = 1.;
    
    hitChi2Max            = _maxchi2dxy;
    
    for (int f=0; f<RobustHelixFinderData::kNTotalFaces; ++f){
      facez         = &helixData._oTracker[f];
      nhitsPerFace  = facez->fNHits;
	  
      if (nhitsPerFace == 0)                       continue;
      bool  isFaceUsed(false);
      RobustHelixFinderData::HitInfo_t  hitInfo;

      for (int i=0; i<nhitsPerFace;++i){
	hit =  &facez->fHitData.at(i);
	if (!hit->_flag.hasAnyProperty(_outlier))   {
	  isFaceUsed = true;
	  break;//skip the faces where there is already a hit
	}
	XYVec rvec = (XYVec(hit->pos().x(),hit->pos().y())-helCenter);
	dr       = sqrtf(rvec.Mag2()) - r;
	wtXY     = _hfit.evalWeightXY(*hit, helCenter);
	drChi2   = sqrtf(dr*dr*wtXY);

	phi_pred = facez->z*dfdz + phi0;
	dphi     = phi_pred - hit->helixPhi();
	wtZPhi   = _hfit.evalWeightZPhi(*hit,helCenter,r);
	dphiChi2 = sqrtf(dphi*dphi*wtZPhi);
	
	hitChi2  = (drChi2 + dphiChi2)/2.;
	
	if ( (drChi2<_maxchi2dxy) && (dphiChi2<_maxchi2dzphi) && (hitChi2 < hitChi2Max)){
	  hitChi2Max  = hitChi2;
	  
	  hitInfo.face          = f;	 
	  hitInfo.faceHitIndex  = i;
	  hitInfo.weightXY      = wtXY;
	  hitInfo.weightZPhi    = wtZPhi;
   
	}
      }//end loop over the hits within the face
      
      if( (!isFaceUsed) && (hitInfo.face >=0 ) ) {
	bestHit.face          = hitInfo.face        ;	 
	bestHit.faceHitIndex  = hitInfo.faceHitIndex;
	bestHit.weightXY      = hitInfo.weightXY    ;
	bestHit.weightZPhi    = hitInfo.weightZPhi  ;
      }
    }//end loop pver the faces
    
    
    if ( (bestHit.face >= 0) ){
      facez  = &helixData._oTracker[bestHit.face];
      hit    = &facez->fHitData.at(bestHit.faceHitIndex);
	

      //add the point 
      hit->_flag.clear(_outlier);

      helixData._sxy.addPoint(hit->pos().x(), hit->pos().y(), bestHit.weightXY);
      helixData._nXYSh   += hit->nStrawHits();

      helixData._szphi.addPoint(facez->z, hit->helixPhi(), bestHit.weightZPhi);
      helixData._nZPhiSh += hit->nStrawHits();


      //update the helix
      updateHelixInfo(helixData);
      
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
    RobustHelixFinderData::FaceZ_t*      facez;

    ComboHit*  hit(0);
    ComboHit*  worsthit(0);

    NRemovedStrawHits = 0;

    while (changed && nhit > 0)
      {
	nhit = 0;
	changed = false;
	accumulator_set<float, stats<tag::median(with_p_square_quantile) > > accx;
	accumulator_set<float, stats<tag::median(with_p_square_quantile) > > accy;

	for (int f=0; f<RobustHelixFinderData::kNTotalFaces; ++f){
	  facez              = &HelixData._oTracker[f];
	  int  nhitsPerFace  = facez->fNHits;
	  
	  if (nhitsPerFace == 0)                     continue;
	  for (int i=0; i<nhitsPerFace;++i){
	    hit =  &facez->fHitData.at(i);
	    bool trashHit=hit->_flag.hasAnyProperty(_outlier);
	    if (trashHit)                              continue;
	    accx(hit->_pos.x());
	    accy(hit->_pos.y());
	    ++nhit;
	  }
	}
	
	float mx = extract_result<tag::median>(accx);
	float my = extract_result<tag::median>(accy);
	float mphi = polyAtan2(my,mx);//atan2f(my,mx);

	float maxdphi{0.0};
	// auto worsthit = hhits.end();
	for (int f=0; f<RobustHelixFinderData::kNTotalFaces; ++f){
	  facez              = &HelixData._oTracker[f];
	  int  nhitsPerFace  = facez->fNHits;
	  
	  if (nhitsPerFace == 0)                     continue;
	  for (int i=0; i<nhitsPerFace;++i){
	    hit =  &facez->fHitData.at(i);
	    bool trashHit = hit->_flag.hasAnyProperty(_outlier);
	    if (trashHit)                              continue;
	    float phi  = polyAtan2(hit->pos().y(), hit->pos().x());//ihit->pos().phi();
	    float dphi = fabs(Angles::deltaPhi(phi,mphi));
	    if(dphi > maxdphi)
	      {
		maxdphi = dphi;
		worsthit = hit;
	      }
	  }
	}

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
    //    const auto& hhits = helixData._hseed.hits();
    accumulator_set<float, stats<tag::weighted_variance(lazy)>, float > terr;
    
    ComboHit*      hit(0);
    RobustHelixFinderData::FaceZ_t*       facez(0);
    int            nhitsFace(0);
    
    for (int f=0; f<RobustHelixFinderData::kNTotalFaces; ++f){
      facez     = &helixData._oTracker[f];
      nhitsFace = facez->fNHits;

      if (nhitsFace == 0)                          continue;

      for (int ip=0; ip<nhitsFace; ++ip){
	hit = &facez->fHitData.at(ip);
	
	if (hit->_flag.hasAnyProperty(_outlier))   continue;
	float wt = std::pow(1.0/_ttcalc.strawHitTimeErr(),2);
	terr(_ttcalc.comboHitTime(*hit),weight=wt);
      }
    }

    if (helixData._hseed.caloCluster().isNonnull())
      {
	float time = _ttcalc.caloClusterTime(*helixData._hseed.caloCluster());
	float wt = std::pow(1.0/_ttcalc.caloClusterTimeErr(helixData._hseed.caloCluster()->diskId()),2);
	terr(time,weight=wt);
      }

    if (sum_of_weights(terr) > 0.0)
      {
	helixData._hseed._t0._t0 = extract_result<tag::weighted_mean>(terr) + _t0shift; // ad-hoc correction FIXME!!
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
    // if (HelixData.shpos() != 0) {
    int loc;
    StrawHitFlag flag;
    for (int i=0; i<size; ++i) {
      loc = shIndices[i];
      const ComboHit ch  = _hfResult._chcol->at(loc);

      if(ch.flag().hasAnyProperty(_hsel) && !ch.flag().hasAnyProperty(_hbkg)  /* && !flag->hasAnyProperty(StrawHitFlag::bkg)*/) {
	ComboHit hhit(ch);
	hhit._flag.clear(StrawHitFlag::resolvedphi);

	cx.Station                 = ch.sid().station();//straw.id().getStation();
	cx.Plane                   = ch.sid().plane() % 2;//straw.id().getPlane() % 2;
	cx.Face                    = ch.sid().face();
	cx.Panel                   = ch.sid().panel();//straw.id().getPanel();

	// get Z-ordered location
	HelixData.orderID(&cx, &co);
     
	int os       = co.Station; 
	int of       = co.Face;
	// int op       = co.Panel;

	int       stationId = os;
	int       faceId    = of + stationId*RobustHelixFinderData::kNFaces;
	// int       panelId   = op + faceId*RobustHelixDataFinderData::kNPanelsPerFace;
	RobustHelixFinderData::FaceZ_t* fz        = &HelixData._oTracker[faceId];

	if ((os < 0) || (os >= RobustHelixFinderData::kNStations     )) printf(" >>> ERROR: wrong station number: %i\n",os);
	if ((of < 0) || (of >= RobustHelixFinderData::kNFaces        )) printf(" >>> ERROR: wrong face    number: %i\n",of);
	//	if ((op < 0) || (op >= RobustHelixFinderData::kNPanelsPerFace)) printf(" >>> ERROR: wrong panel   number: %i\n",op);

	// pz->fHitData.push_back(RobustHelixDataPoint(loc,sh,shp,straw,flag));
	fz->fHitData.push_back(hhit);
	fz->fNHits  = fz->fNHits + 1;
	if (fz->fNHits > RobustHelixFinderData::kNMaxHitsPerFace) printf("[RobustHelixDataFinderAlg::fillFaceOrderedHits] number of hits with the panel exceed the limit: NHits =  %i MaxNHits = %i\n", fz->fNHits, RobustHelixFinderData::kNMaxHitsPerFace);
	++nFiltComboHits;
	nFiltStrawHits += ch.nStrawHits();
      }
    }
    // }
    
    HelixData._nFiltComboHits = nFiltComboHits;  //ComboHit counter
    HelixData._nFiltStrawHits = nFiltStrawHits;  //StrawHit counter
    
  }

  unsigned  RobustHelixFinder::filterCircleHits(RobustHelixFinderData& helixData)
  {
    unsigned changed(0);
    int      nGoodSH(0);
    static XYZVec zaxis(0.0,0.0,1.0); // unit in z direction
    RobustHelix& helix = helixData._hseed._helix;

    // loop over hits
    ComboHit*     hit(0);
    RobustHelixFinderData::FaceZ_t*      facez;
    int           nhitsFace(0);

    for (int f=0; f<RobustHelixFinderData::kNTotalFaces; ++f){
      facez     = &helixData._oTracker[f];
      nhitsFace = facez->fNHits;
      if (nhitsFace == 0)                        continue;

      for (int ip=0; ip<nhitsFace; ++ip){
	hit = &facez->fHitData.at(ip);
	
	bool oldout = hit->_flag.hasAnyProperty(_outlier);
	hit->_flag.clear(_outlier);

	const XYZVec& wdir = hit->wdir();
	XYZVec cvec = PerpVector(hit->pos() - helix.center(),Geom::ZDir()); // direction from the circle center to the hit
	XYZVec cdir = cvec.Unit(); // direction from the circle center to the hit
	float rwdot = wdir.Dot(cdir); // compare directions of radius and wire
	if(rwdot > _maxrwdot){
	  hit->_flag.merge(_outlier);
	  if(!oldout) ++changed;
	  continue;
	}
	float dr = sqrtf(cvec.mag2())-helix.radius();
	if ( fabs(dr) > _maxdr ) {
	  hit->_flag.merge(_outlier);
	  if(!oldout) ++changed;
	  continue;
	}
	
	float rwdot2 = rwdot*rwdot;
	// compute radial difference and pull
	float werr = hit->posRes(StrawHitPosition::wire);
	float terr = hit->posRes(StrawHitPosition::trans);
	// the resolution is dominated the resolution along the wire
	float rres = std::max(sqrtf(werr*werr*rwdot2 + terr*terr*(1.0-rwdot2)),_minrerr);
	float rpull = dr/rres;
	if ( fabs(rpull) > _maxrpull ) {
	  hit->_flag.merge(_outlier);
	  if(!oldout) ++changed;
	  continue;
	}
	
	if (oldout) ++changed;
	
	nGoodSH += hit->nStrawHits();
      }//end loop over the hits within a face
    }//end loop over the faces
    
    helixData._nXYSh = nGoodSH;
    
    return changed;
  }


  void RobustHelixFinder::fitHelix(RobustHelixFinderData& helixData){
    // iteratively fit the helix including filtering
    unsigned niter(0);
    unsigned nitermva(0);
    bool changed(true);

    do {
      unsigned xyniter(0);
      unsigned xychanged = filterCircleHits(helixData);
      while (helixData._hseed._status.hasAllProperties(TrkFitFlag::circleOK) && xyniter < _maxniter && xychanged) {
	_hfit.fitCircle(helixData);
	xychanged = filterCircleHits(helixData);
	++xyniter;
      } 
   
      if (_diag) helixData._diag.xyniter = xyniter;

      // then fit phi-Z
      if (helixData._hseed._status.hasAnyProperty(TrkFitFlag::circleOK)) {
	if (xyniter < _maxniter)
	  helixData._hseed._status.merge(TrkFitFlag::circleConverged);
	else
	  helixData._hseed._status.clear(TrkFitFlag::circleConverged);

	// solve for the longitudinal parameters
	unsigned fzniter(0);
	bool fzchanged(false);
	do {
	  _hfit.fitFZ(helixData);
	  fzchanged = filterHits(helixData);
	  ++fzniter;
	} while (helixData._hseed._status.hasAllProperties(TrkFitFlag::phizOK)  && fzniter < _maxniter && fzchanged);

	if (_diag) helixData._diag.fzniter = fzniter;

	if (helixData._hseed._status.hasAnyProperty(TrkFitFlag::phizOK)) {
	  if (fzniter < _maxniter)
	    helixData._hseed._status.merge(TrkFitFlag::phizConverged);
	  else
	    helixData._hseed._status.clear(TrkFitFlag::phizConverged);
	}
      }
      //here is where we should check for the hits within the face to searchfor missing/best ones
      ++niter;

      // update the stereo hit positions; this checks how much the positions changed
      // do this only in non trigger mode

      if (_updateStereo && _hfit.goodHelix(helixData._hseed.helix())) changed = updateStereo(helixData);      
    } while (_hfit.goodHelix(helixData._hseed.helix()) && niter < _maxniter && changed);

    if (_diag) helixData._diag.niter = niter;

    if (_hfit.goodHelix(helixData._hseed.helix())) {
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
      _nitermva->Fill(nitermva);
      if (!_usemva) fillMVA(helixData);
    }
  }



  //------------------------------------------------------------------------------------------


  void RobustHelixFinder::fitHelix_2(RobustHelixFinderData& helixData){
    // iteratively fit the helix including filtering
    unsigned nitermva(0);

    unsigned xyniter(0);
    int      xychanged = filterXYHits(helixData);
    while (helixData._hseed._status.hasAnyProperty(TrkFitFlag::circleOK) && xyniter < _maxniter && (xychanged!=0)) {
      //	_hfit.fitCircle(helixData);//now the filterCircleHits updates also the info of the helixs
      xychanged = filterXYHits(helixData);
      ++xyniter;
    } 
   
    if (_diag) helixData._diag.xyniter = xyniter;

    // then fit phi-Z
    if (helixData._hseed._status.hasAnyProperty(TrkFitFlag::circleOK)) {
      if (xyniter < _maxniter)
	helixData._hseed._status.merge(TrkFitFlag::circleConverged);
      else
	helixData._hseed._status.clear(TrkFitFlag::circleConverged);

      // solve for the longitudinal parameters
      unsigned fzniter(0);
      _hfit.fitFZ_2(helixData);
      int fzchanged = filterZPhiHits(helixData);//(false);
      while (helixData._hseed._status.hasAnyProperty(TrkFitFlag::phizOK)  && fzniter < _maxniter && (fzchanged!=0)) {
	fzchanged = filterZPhiHits(helixData);
	++fzniter;
      } //while (helixData._hseed._status.hasAllProperties(TrkFitFlag::phizOK)  && fzniter < _maxniter && fzchanged);

      if (_diag) helixData._diag.fzniter = fzniter;

      if (helixData._hseed._status.hasAnyProperty(TrkFitFlag::phizOK)) {
	if (fzniter < _maxniter)
	  helixData._hseed._status.merge(TrkFitFlag::phizConverged);
	else
	  helixData._hseed._status.clear(TrkFitFlag::phizConverged);
      }
    }

    //now search for missing hits
    findMissingHits(helixData);

    if (_hfit.goodHelix(helixData._hseed.helix()) && _hfit.goodHelixChi2(helixData)) {
      helixData._hseed._status.merge(TrkFitFlag::helixOK);
      updateT0(helixData);
      
      helixData._hseed._status.merge(TrkFitFlag::helixConverged);

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
      //      _niter->Fill(niter);
      _nitermva->Fill(nitermva);
      if (!_usemva) fillMVA(helixData);
    }
  }

  void RobustHelixFinder::refitHelix(RobustHelixFinderData& helixData) {
    // reset the fit status flags, in case this is called iteratively
    helixData._hseed._status.clear(TrkFitFlag::helixOK);      
    _hfit.fitCircle(helixData);
    if (helixData._hseed._status.hasAnyProperty(TrkFitFlag::circleOK)) {
      _hfit.fitFZ(helixData);
      if (_hfit.goodHelix(helixData._hseed._helix)) helixData._hseed._status.merge(TrkFitFlag::helixOK);
    }
  }

  unsigned RobustHelixFinder::hitCount(RobustHelixFinderData& helixData) {
    // return std::count_if(helixData._hseed._hhits.begin(),helixData._hseed._hhits.end(),
    // 			 [&](const ComboHit& hhit){return !hhit.flag().hasAnyProperty(_outlier);});
    unsigned nHits(0);
    
    ComboHit*      hit(0);
    RobustHelixFinderData::FaceZ_t*       facez(0);
    int            nhitsFace(0);
    
    for (int f=0; f<RobustHelixFinderData::kNTotalFaces; ++f){
      facez     = &helixData._oTracker[f];
      nhitsFace = facez->fNHits;

      if (nhitsFace == 0)                          continue;

      for (int ip=0; ip<nhitsFace; ++ip){
	hit = &facez->fHitData.at(ip);
	
	if (hit->_flag.hasAnyProperty(_outlier))   continue;

	++nHits;
      }
    }

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
    double          mm2MeV = 3./10.;  // approximately , at B=1T
	  
    int loc = _data.nseeds[helCounter];
    if (loc < _data.maxSeeds()) {
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
   
      _data.nfz0counter [helCounter][loc] = helixData._diag.nfz0counter;

      _data.nshsxy_0    [helCounter][loc] = helixData._diag.nshsxy_0;
      _data.rsxy_0      [helCounter][loc] = helixData._diag.rsxy_0;
      _data.chi2dsxy_0  [helCounter][loc] = helixData._diag.chi2dsxy_0;
   	                
      _data.nshsxy_1    [helCounter][loc] = helixData._diag.nshsxy_1;
      _data.rsxy_1      [helCounter][loc] = helixData._diag.rsxy_1;
      _data.chi2dsxy_1  [helCounter][loc] = helixData._diag.chi2dsxy_1;
   	                
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
	
      // _data.chi2XY[loc]   = helixData._sxy.chi2DofCircle();
      // _data.chi2ZPhi[loc] = helixData._szphi.chi2DofLine();
	    
      _data.nseeds[helCounter]++;
	    
      _data.dr           [helCounter][loc] = helixData._diag.radius_2 - helixData._diag.radius_1;
      _data.chi2d_helix  [helCounter][loc] = helixData._diag.chi2d_helix;
    }   else {
      printf(" N(seeds) > %i, IGNORE SEED\n",_data.maxSeeds());
    }
  }

  void     RobustHelixFinder::updateHelixInfo(RobustHelixFinderData& helixData){
    RobustHelix&  helix        = helixData._hseed._helix;
    XYVec         center       = XYVec(helix.center().x(), helix.center().y());
    float         radius       = helixData._sxy.radius();
    center.SetX(helixData._sxy.x0());
    center.SetY(helixData._sxy.y0());

    //update the LSqsums
    ComboHit*      hit(0);
    RobustHelixFinderData::FaceZ_t*       facez(0);
    int            nhitsFace(0);
    float          wtXY, wtZPhi;
    
    helixData._sxy.clear();
    helixData._szphi.clear();
    helixData._nXYSh   = 0;
    helixData._nZPhiSh = 0;

    for (int f=0; f<RobustHelixFinderData::kNTotalFaces; ++f){
      facez     = &helixData._oTracker[f];
      nhitsFace = facez->fNHits;

      if (nhitsFace == 0)                          continue;

      for (int ip=0; ip<nhitsFace; ++ip){
	hit    = &facez->fHitData.at(ip);
	 
	if (hit->_flag.hasAnyProperty(_outlier))   continue;
	
	wtXY   = _hfit.evalWeightXY(*hit, center);
	wtZPhi = _hfit.evalWeightZPhi(*hit, center, radius);

	helixData._sxy.addPoint(hit->pos().x(), hit->pos().y(), wtXY);
	helixData._szphi.addPoint(facez->z, hit->helixPhi(), wtZPhi);

	helixData._nXYSh   += hit->nStrawHits();
	helixData._nZPhiSh += hit->nStrawHits();
      }
      //should we update the parameters (center, radius) after a while?
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
  void     RobustHelixFinder::searchWorstHitXY(RobustHelixFinderData& helixData, RobustHelixFinderData::HitInfo_t& hitInfo){
    hitInfo.face         = -1;
    hitInfo.faceHitIndex = -1;
    
    float          dr, wt, hitChi2;
    XYVec          centerXY     = XYVec(helixData._hseed._helix.center().x(), helixData._hseed._helix.center().y());
    XYZVec         center       = helixData._hseed._helix.center();
    float          helix_radius = helixData._hseed._helix.radius();
    float          hitChi2Worst = _maxrpull*_maxrpull;

    ComboHit*      hit(0);
    RobustHelixFinderData::FaceZ_t*       facez(0);
    int            nhitsFace(0);
    
    for (int f=0; f<RobustHelixFinderData::kNTotalFaces; ++f){
      facez     = &helixData._oTracker[f];
      nhitsFace = facez->fNHits;

      if (nhitsFace == 0)                          continue;

      for (int ip=0; ip<nhitsFace; ++ip){
	hit = &facez->fHitData.at(ip);
	
	if (hit->_flag.hasAnyProperty(_outlier))   continue;
	
	XYZVec  cvec  = PerpVector(hit->pos() - center,Geom::ZDir());
	dr    = sqrtf(cvec.mag2()) - helix_radius;
	wt    = _hfit.evalWeightXY(*hit, centerXY);

	hitChi2 = dr*dr*wt;

	if (hitChi2 > hitChi2Worst) {
	  hitChi2Worst         = hitChi2;
	  hitInfo.face         = f;
	  hitInfo.faceHitIndex = ip;
	  hitInfo.weightXY     = wt;
	}
      }
    }
  }
  
  //----------------------------------------------------------------------------------------------------
  void     RobustHelixFinder::searchWorstHitZPhi(RobustHelixFinderData& helixData, RobustHelixFinderData::HitInfo_t& hitInfo){
    hitInfo.face         = -1;
    hitInfo.faceHitIndex = -1;
    
    // XYZVec         center       = helixData._hseed._helix.center();
    XYVec          centerXY     = XYVec(helixData._hseed._helix.center().x(), helixData._hseed._helix.center().y());
    float          helix_radius = helixData._hseed._helix.radius();

    ComboHit*      hit(0);
    RobustHelixFinderData::FaceZ_t*       facez(0);
    int            nhitsFace(0);
    
    float          chi2min(1e10), chi2;
    ::LsqSums4     szphi;

    for (int f=0; f<RobustHelixFinderData::kNTotalFaces; ++f){
      facez     = &helixData._oTracker[f];
      nhitsFace = facez->fNHits;

      if (nhitsFace == 0)                          continue;

      for (int ip=0; ip<nhitsFace; ++ip){
	hit = &facez->fHitData.at(ip);
	
	if (hit->_flag.hasAnyProperty(_outlier))   continue;
	
	szphi.init(helixData._szphi);

	// XYZVec  cvec  = PerpVector(hit->pos() - center,Geom::ZDir());
	float   phi   = hit->helixPhi();
	float   wt    = _hfit.evalWeightZPhi(*hit, centerXY, helix_radius);

	szphi.removePoint(facez->z, phi, wt);
	chi2 = szphi.chi2DofLine();
	

	if (chi2 < chi2min) {
	  chi2min              = chi2;
	  hitInfo.face         = f;
	  hitInfo.faceHitIndex = ip;
	  hitInfo.weightZPhi   = wt;
	}
      }
    }
  }
}
using mu2e::RobustHelixFinder;
DEFINE_ART_MODULE(RobustHelixFinder);
