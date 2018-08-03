///////////////////////////////////////////////////////////////////////////////
// Calorimeter-driven track finding
// Pattern recognition only, passes results to CalSeedFit
// P.Murat, G.Pezzullo
// try to order routines alphabetically
///////////////////////////////////////////////////////////////////////////////
#include "fhiclcpp/ParameterSet.h"

#include "CalPatRec/inc/CalHelixFinder_module.hh"

// framework
#include "art/Framework/Principal/Handle.h"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"

// conditions
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
// #include "CalPatRec/inc/KalFitResult.hh"
#include "RecoDataProducts/inc/StrawHitIndex.hh"

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/algorithm/string.hpp>

#include "CalPatRec/inc/CalHelixFinderData.hh"

#include "Mu2eUtilities/inc/ModuleHistToolBase.hh"
#include "art/Utilities/make_tool.h"
#include "Mu2eUtilities/inc/polyAtan2.hh"

#include "TVector2.h"
#include "TSystem.h"
#include "TInterpreter.h"

using namespace std;
using namespace boost::accumulators;
using CLHEP::HepVector;
using CLHEP::Hep3Vector;

namespace mu2e {
//-----------------------------------------------------------------------------
// module constructor, parameter defaults are defiend in CalPatRec/fcl/prolog.fcl
//-----------------------------------------------------------------------------
  CalHelixFinder::CalHelixFinder(fhicl::ParameterSet const& pset) :
    _diagLevel   (pset.get<int>   ("diagLevel"                      )),
    _debugLevel  (pset.get<int>   ("debugLevel"                     )),
    _printfreq   (pset.get<int>   ("printFrequency"                 )),
    _useAsFilter (pset.get<int>   ("useAsFilter"                    )),
    _shLabel     (pset.get<string>("StrawHitCollectionLabel"        )),
    // _shpLabel    (pset.get<string>("StrawHitPositionCollectionLabel")),
    _shfLabel    (pset.get<string>("StrawHitFlagCollectionLabel"    )),
    _timeclLabel (pset.get<string>("TimeClusterCollectionLabel"       )),
    _minNHitsTimeCluster(pset.get<int>("minNHitsTimeCluster"       )),
    _tpart       ((TrkParticle::type)(pset.get<int>("fitparticle"))),
    _fdir        ((TrkFitDirection::FitDirection)(pset.get<int>("fitdirection"))),
    _hfinder     (pset.get<fhicl::ParameterSet>("HelixFinderAlg",fhicl::ParameterSet()))
  {
    consumes<ComboHitCollection>(_shLabel);
    consumes<StrawHitFlagCollection>(_shfLabel);
    consumes<TimeClusterCollection>(_timeclLabel);
    produces<HelixSeedCollection>();
//-----------------------------------------------------------------------------
// provide for interactive disanostics
//-----------------------------------------------------------------------------
    _helTraj          = 0;
    _timeOffsets      = new fhicl::ParameterSet(pset.get<fhicl::ParameterSet>("TimeOffsets",fhicl::ParameterSet()));

    _data.shLabel     = _shLabel;
    _data.timeOffsets = _timeOffsets;

    if (_debugLevel != 0) _printfreq = 1;

    if (_diagLevel != 0) _hmanager = art::make_tool<ModuleHistToolBase>(pset.get<fhicl::ParameterSet>("diagPlugin"));
    else                 _hmanager = std::make_unique<ModuleHistToolBase>();

  }

//-----------------------------------------------------------------------------
// destructor
//-----------------------------------------------------------------------------
  CalHelixFinder::~CalHelixFinder() {
    if (_helTraj) delete _helTraj;
    delete _timeOffsets;
  }

//-----------------------------------------------------------------------------
  void CalHelixFinder::beginJob(){
    art::ServiceHandle<art::TFileService> tfs;
    _hmanager->bookHistograms(tfs);
  }

//-----------------------------------------------------------------------------
  bool CalHelixFinder::beginRun(art::Run& ) {
    mu2e::GeomHandle<mu2e::TTracker> th;
    _tracker = th.get();

    mu2e::GeomHandle<mu2e::Calorimeter> ch;
    _calorimeter = ch.get();

    _hfinder.setTracker    (_tracker);
    _hfinder.setCalorimeter(_calorimeter);

    ChannelID cx, co;
    // int       nTotalStations = _tracker->nStations();

    for (int ist=0; ist<_tracker->nStations(); ist++) {
      const Station* st = &_tracker->getStation(ist);

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
	  int op = co.Panel;

	  int       stationId = os;
	  int       faceId    = of + stationId*StrawId::_nfaces*FaceZ_t::kNPlanesPerStation;//CalHelixFinderData::kNFaces;
	  //	  int       panelId   = op;// + faceId*CalHelixFinderData::kNPanelsPerFace;
	  FaceZ_t*  fz        = &_hfResult._oTracker[faceId];
	  fz->z               = (panel->getStraw(0).getMidPoint().z()+panel->getStraw(1).getMidPoint().z())/2.;

	  PanelZ_t* pz        = &fz->panelZs[op];//_hfResult._oTracker[panelId];

	  //	  pz->fPanel = panel;
	  //-----------------------------------------------------------------------------
	  // panel caches phi of its center and the z
	  //-----------------------------------------------------------------------------
	  pz->wx     = panel->straw0Direction().x();
	  pz->wy     = panel->straw0Direction().y();
	  pz->phi    = TVector2::Phi_0_2pi(polyAtan2(panel->straw0MidPoint().y(),panel->straw0MidPoint().x()));
	  pz->fNHits = 0;
	}	
      }
    }

    if (_debugLevel > 10){
      printf("//-----------------------------------//\n");
      printf("//     Face      Panel      Z        //\n");
      printf("//-----------------------------------//\n");

      FaceZ_t*        facez(0);
      //      PanelZ_t*       panelz(0);

      for (int f=0; f<StrawId::_ntotalfaces; ++f){
	facez     = &_hfResult._oTracker[f];
	double z  = facez->z;
	for (int p=0; p<FaceZ_t::kNPanels; ++p){
	  //	  panelz =  &facez->panelZs[p];//	panelz = &_hfResult._oTracker[p];
	  printf("//  %5i      %5i     %10.3f //\n", f, p, z);
	}
      }
      printf("//----------------------------------//\n");

    }

    return true;
  }

//-----------------------------------------------------------------------------
// find the input data objects
//-----------------------------------------------------------------------------
  bool CalHelixFinder::findData(const art::Event& evt) {

    if (evt.getByLabel(_shLabel, _strawhitsH)) {
      _chcol = _strawhitsH.product();
    }
    else {
      _chcol  = 0;
      printf(" >>> ERROR in CalHelixFinder::findData: StrawHitCollection with label=%s not found.\n",
             _shLabel.data());
    }

    // art::Handle<mu2e::StrawHitPositionCollection> shposH;
    // if (evt.getByLabel(_shpLabel,shposH)) {
    //   _shpcol = shposH.product();
    // }
    // else {
    //   _shpcol = 0;
    //   printf(" >>> ERROR in CalHelixFinder::findData: StrawHitPositionCollection with label=%s not found.\n",
    //          _shpLabel.data());
    // }

    art::Handle<mu2e::StrawHitFlagCollection> shflagH;
    if (evt.getByLabel(_shfLabel,shflagH)) {
      _shfcol = shflagH.product();
    }
    else {
      _shfcol = 0;
      printf(" >>> ERROR in CalHelixFinder::findData: StrawHitFlagCollection with label=%s not found.\n",
             _shfLabel.data());
    }


    if (evt.getByLabel(_timeclLabel, _timeclcolH)) {
      _timeclcol = _timeclcolH.product();
    }
    else {
      _timeclcol = 0;
      printf(" >>> ERROR in CalHelixFinder::findData: TimeClusterCollection with label=%s not found.\n",
             _timeclLabel.data());
    }
//-----------------------------------------------------------------------------
// done
//-----------------------------------------------------------------------------
    return (_chcol != 0) && (_shfcol != 0) /*&& (_shpcol != 0) */&& (_timeclcol != 0);
  }

//-----------------------------------------------------------------------------
// event entry point
//-----------------------------------------------------------------------------
  bool CalHelixFinder::filter(art::Event& event ) {
    const char*             oname = "CalHelixFinder::filter";
    //    CalHelixFinderData      hf_result;
                                        // diagnostic info
    _data.event     = &event;
    _data.nseeds[0] = 0;
    _data.nseeds[1] = 0;
    _iev            = event.id().event();
    int   nGoodTClusterHits(0);

    if ((_debugLevel > 0) && (_iev%_printfreq) == 0) printf("[%s] : START event number %8i\n", oname,_iev);

    unique_ptr<HelixSeedCollection>    outseeds(new HelixSeedCollection);
//-----------------------------------------------------------------------------
// find the data
//-----------------------------------------------------------------------------
    if (!findData(event)) {
      printf("%s ERROR: No straw hits found, RETURN\n", oname);
                                                            goto END;
    }
//-----------------------------------------------------------------------------
// loop over found time peaks - for us, - "eligible" calorimeter clusters
//-----------------------------------------------------------------------------
    _hfResult._tpart  = _tpart;
    _hfResult._fdir   = _fdir;
    _hfResult._chcol  = _chcol;
    // _hfResult._shpos  = _shpcol;
    _hfResult._shfcol = _shfcol;

    _data.nTimePeaks  = _timeclcol->size();
    for (int ipeak=0; ipeak<_data.nTimePeaks; ipeak++) {
      const TimeCluster* tc = &_timeclcol->at(ipeak);
      nGoodTClusterHits     = goodHitsTimeCluster(tc);
      if ( nGoodTClusterHits < _minNHitsTimeCluster)         continue;

      HelixSeed          helix_seed;
//-----------------------------------------------------------------------------
// create track definitions for the helix fit from this initial information
// track fitting objects for this peak
//-----------------------------------------------------------------------------
      _hfResult.clearTempVariables();

      _hfResult._timeCluster    = tc;
      _hfResult._timeClusterPtr = art::Ptr<mu2e::TimeCluster>(_timeclcolH,ipeak);
//-----------------------------------------------------------------------------
// Step 1: pattern recognition. Find initial helical approximation of a track
//-----------------------------------------------------------------------------
      int rc = _hfinder.findHelix(_hfResult,tc);

      if (rc) {
//-----------------------------------------------------------------------------
// fill seed information
//-----------------------------------------------------------------------------
        initHelixSeed(helix_seed, _hfResult);
        helix_seed._status.merge(TrkFitFlag::helixOK);
        outseeds->push_back(helix_seed);

        if (_diagLevel > 0) {
//--------------------------------------------------------------------------------
// fill diagnostic information
//--------------------------------------------------------------------------------
          int             nhitsMin(15);
          double          mm2MeV = 3./10.;  // approximately , at B=1T

          int loc = _data.nseeds[0];
          if (loc < _data.maxSeeds()) {
            int nhits          = helix_seed._hhits.size();
            _data.ntclhits[loc]= nGoodTClusterHits;
            _data.nhits[loc]   = nhits;
            _data.radius[loc]  = helix_seed.helix().radius();
            _data.pT[loc]      = mm2MeV*_data.radius[loc];
            _data.p[loc]       = _data.pT[loc]/std::cos( std::atan(helix_seed.helix().lambda()/_data.radius[loc]));

            _data.chi2XY[loc]   = _hfResult._sxy.chi2DofCircle();
            _data.chi2ZPhi[loc] = _hfResult._szphi.chi2DofLine();

            _data.nseeds[0]++;
            _data.good[loc] = 0;
            if (nhits >= nhitsMin) {
              _data.nseeds[1]++;
              _data.good[loc] = 1;
            }
            _data.nStationPairs[loc] = _hfResult._diag.nStationPairs;

            _data.dr           [loc] = _hfResult._diag.dr;
            _data.shmeanr      [loc] = _hfResult._diag.straw_mean_radius;
            _data.chi2d_helix  [loc] = _hfResult._diag.chi2d_helix;
            if (_hfResult._diag.chi2d_helix>3) printf("[%s] : chi2Helix = %10.3f event number %8i\n", oname,_hfResult._diag.chi2d_helix,_iev);
//-----------------------------------------------------------------------------
// info of the track candidate after the first loop with findtrack on CalHelixFinderAlg::doPatternRecognition
//-----------------------------------------------------------------------------
            _data.loopId       [loc] = _hfResult._diag.loopId_4;
            if (_hfResult._diag.loopId_4 == 1) {
              _data.chi2d_loop0       [loc] = _hfResult._diag.chi2_dof_circle_12;
              _data.chi2d_line_loop0  [loc] = _hfResult._diag.chi2_dof_line_13;
              _data.npoints_loop0     [loc] = _hfResult._diag.n_active_11;

            }
            if (_hfResult._diag.loopId_4 == 2){
              _data.chi2d_loop1       [loc] = _hfResult._diag.chi2_dof_circle_12;
              _data.chi2d_line_loop1  [loc] = _hfResult._diag.chi2_dof_line_13;
              _data.npoints_loop1     [loc] = _hfResult._diag.n_active_11;
            }

//--------------------------------------------------------------------------------
// info of the track candidate during the CAlHelixFinderAlg::findTrack loop
//--------------------------------------------------------------------------------
	    int   counter(0);
	    for (int f=0; f<StrawId::_ntotalfaces; ++f){
	      FaceZ_t* facez     = &_hfResult._oTracker[f];
	      for (int p=0; p<FaceZ_t::kNPanels; ++p){//for (int p=0; p<CalHelixFinderData::kNTotalPanels; ++p){
		PanelZ_t* panelz = &facez->panelZs[p];//&_hfResult._oTracker[p];
		int       nhits  = panelz->fNHits;
		if (nhits == 0)                                  continue;
	      
		for (int i=0; i<nhits; ++i){   
		  //		  ComboHit*	hit = &panelz->fHitData.at(i);
		  int index = facez->evalUniqueHitIndex(f,p,i);//p*CalHelixFinderData::kNMaxHitsPerPanel + i;
		  if (_hfResult._hitsUsed[index] != 1)           continue;
		
		  // double   dzFromSeed = hit->_dzFromSeed;     //distance form the hit used to seed the 3D-search
		  // double   drFromPred = hit->_drFromPred;     //distance from prediction
		  // _data.hitDzSeed[loc][counter] = dzFromSeed;
		  // _data.hitDrPred[loc][counter] = drFromPred;
		  ++counter;
		}//end loop over the hits within a panel
	      }//end panels loop
	    }//end faces loop
	  }
	  else {
	    printf(" N(seeds) > %i, IGNORE SEED\n",_data.maxSeeds());
	  }
	}
      }
    }
//--------------------------------------------------------------------------------
// fill histograms
//--------------------------------------------------------------------------------
    if (_diagLevel > 0) _hmanager->fillHistograms(&_data);
//-----------------------------------------------------------------------------
// put reconstructed tracks into the event record
//-----------------------------------------------------------------------------
  END:;
    int    nseeds = outseeds->size();
    event.put(std::move(outseeds));
//-----------------------------------------------------------------------------
// filtering
//-----------------------------------------------------------------------------
    if (_useAsFilter == 0) return true;
    else                   return (nseeds >  0);
 }

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  void CalHelixFinder::endJob() {
    // does this cause the file to close?
    art::ServiceHandle<art::TFileService> tfs;
  }

//--------------------------------------------------------------------------------
// set helix parameters
//-----------------------------------------------------------------------------
  void CalHelixFinder::initHelixSeed(HelixSeed& HelSeed, CalHelixFinderData& HfResult) {

    HelixTraj* hel = HfResult.helix();

    double   helixRadius   = 1./fabs(hel->omega());
    double   impactParam   = hel->d0();
    double   phi0          = hel->phi0();
    double   x0            = -(helixRadius + impactParam)*sin(phi0);
    double   y0            =  (helixRadius + impactParam)*cos(phi0);
    double   dfdz          = 1./hel->tanDip()/helixRadius;
    double   z0            = hel->z0();
                                        // center of the helix in the transverse plane
    Hep3Vector center(x0, y0, 0);
                                        //define the reconstructed helix parameters

    HelSeed._helix._rcent    = center.perp();
    HelSeed._helix._fcent    = center.phi();
    HelSeed._helix._radius   = helixRadius;
    HelSeed._helix._lambda   = 1./dfdz;
    HelSeed._helix._fz0      = -z0*dfdz*_hfinder._dfdzsign + phi0 - M_PI/2.;
    HelSeed._helix._helicity = _hfinder._dfdzsign > 0 ? Helicity::poshel : Helicity::neghel;

                                        //now evaluate the helix T0 using the calorimeter cluster
    double   mm2MeV        = 3/10.;//FIX ME!
    double   tandip        = hel->tanDip();
    double   mom           = helixRadius*mm2MeV/std::cos( std::atan(tandip));
    double   beta          = _tpart.beta(mom);
    CLHEP::Hep3Vector        gpos = _calorimeter->geomUtil().diskToMu2e(HfResult._timeClusterPtr->caloCluster()->diskId(),
                                                                        HfResult._timeClusterPtr->caloCluster()->cog3Vector());
    CLHEP::Hep3Vector        tpos = _calorimeter->geomUtil().mu2eToTracker(gpos);
    double   pitchAngle    = M_PI/2. - atan(tandip);
    double   hel_t0        = HfResult._timeClusterPtr->caloCluster()->time() - (tpos.z() - z0)/sin(pitchAngle)/(beta*CLHEP::c_light);

    HelSeed._t0            = TrkT0(hel_t0, 0.1); //dummy error on T0 FIXME!
    HelSeed._timeCluster   = HfResult._timeClusterPtr;

    // cluster hits assigned to the reconsturcted Helix

    int nhits = HfResult.nGoodHits();
    // printf("[CalHelixFinder::initHelixSeed] radius = %2.3f x0 = %2.3f y0 = %2.3f dfdz = %2.3e nhits = %i chi2XY = %2.3f chi2PHIZ = %2.3f\n",
    //     helixRadius, center.x(), center.y(), dfdz, nhits, HfResult._sxyw.chi2DofCircle(), HfResult._srphi.chi2DofLine());
    // printf("[CalHelixFinder::initHelixSeed] Index      X          Y         Z          PHI\n");

    double     z_start(0);
    HelSeed._hhits.setParent(_chcol->parent());
    for (int i=0; i<nhits; ++i){
      HitInfo_t*      hitInfo = &HfResult._goodhits[i];
      FaceZ_t*        facez   = &HfResult._oTracker[hitInfo->face];
      PanelZ_t*       panelz  = &facez->panelZs[hitInfo->panel];
      
      ComboHit*       hit    = &panelz->fHitData.at(hitInfo->panelHitIndex);

      double                  hit_z  = hit->pos().z();
      if ( i==0 ) z_start = hit_z;

      double                  shphi  = XYZVec(hit->pos() - HelSeed._helix.center()).phi();
      int                     nLoops = (hit_z - z_start)/(2.*M_PI/dfdz);
      shphi = shphi + double(nLoops)*2.*M_PI;

      ComboHit                hhit(*hit);
      hhit._hphi = shphi;
      hhit._flag.merge(StrawHitFlag::resolvedphi);

      HelSeed._hhits.push_back(hhit);
    }
  }

//-----------------------------------------------------------------------------
  int CalHelixFinder::initHelixFinderData(CalHelixFinderData&                Data,
                                          const TrkParticle&                 TPart,
                                          const TrkFitDirection&             FDir,
                                          const ComboHitCollection*          ComboCollection ,
                                          // const StrawHitPositionCollection*  ShPosCollection ,
                                          const StrawHitFlagCollection*      ShFlagCollection) {
    Data._fit         = TrkErrCode::fail;
    Data._tpart       = TPart;
    Data._fdir        = FDir;

    Data._chcol       = ComboCollection;
    // Data._shpos       = ShPosCollection;
    Data._shfcol      = ShFlagCollection;

    Data._radius      = -1.0;
    Data._dfdz        = 0.;
    Data._fz0         = 0.;

    return 0;
  }

  int  CalHelixFinder::goodHitsTimeCluster(const TimeCluster* TCluster){
    int   nhits         = TCluster->nhits();
    int   ngoodhits(0);
    //    std::vector<string> bkgsel;
    //    bkgsel.push_back("Background");
    double     minT(500.), maxT(2000.);
    for (int i=0; i<nhits; ++i){
      int          index   = TCluster->hits().at(i);
      StrawHitFlag flag    = _shfcol->at(index);
      ComboHit     sh      = _chcol ->at(index);
      int          bkg_hit = flag.hasAnyProperty(StrawHitFlag::bkg);
      if (bkg_hit)                              continue;
      if ( (sh.time() < minT) || (sh.time() > maxT) )  continue;

      // ++ngoodhits;
      ngoodhits += sh.nStrawHits();
    }

    return ngoodhits;
  }



}

using mu2e::CalHelixFinder;
DEFINE_ART_MODULE(CalHelixFinder);
