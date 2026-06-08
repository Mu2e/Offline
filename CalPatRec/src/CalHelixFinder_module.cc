///////////////////////////////////////////////////////////////////////////////
// Calorimeter-driven track finding
// Pattern recognition only, passes results to CalSeedFit
// P.Murat, G.Pezzullo
// try to order routines alphabetically
///////////////////////////////////////////////////////////////////////////////
#include "fhiclcpp/ParameterSet.h"

#include "Offline/CalPatRec/inc/CalHelixFinder_module.hh"

// framework
#include "art/Framework/Principal/Handle.h"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "art_root_io/TFileService.h"

// conditions
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/BFieldGeom/inc/BFieldManager.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Offline/RecoDataProducts/inc/StrawHitIndex.hh"
#include "Offline/RecoDataProducts/inc/HelixHit.hh"

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/algorithm/string.hpp>

#include "Offline/CalPatRec/inc/CalHelixFinderData.hh"

#include "Offline/Mu2eUtilities/inc/ModuleHistToolBase.hh"
#include "art/Utilities/make_tool.h"
#include "Offline/Mu2eUtilities/inc/polyAtan2.hh"
#include "Offline/Mu2eUtilities/inc/HelixTool.hh"

#include "TVector2.h"
#include "TSystem.h"
#include "TInterpreter.h"

#include <unordered_set>

using namespace boost::accumulators;
using CLHEP::HepVector;
using CLHEP::Hep3Vector;

namespace mu2e {
  //-----------------------------------------------------------------------------
  // module constructor, parameter defaults are defiend in CalPatRec/fcl/prolog.fcl
  //-----------------------------------------------------------------------------
  CalHelixFinder::CalHelixFinder(const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config},
    _diagLevel(config().diagLevel()),
    _debugLevel(config().debugLevel()),
    _printfreq(config().printfreq()),
    _shLabel(config().shLabel()),
    _timeclLabel(config().timeclLabel()),
    _minNHitsTimeCluster(config().minNHitsTimeCluster()),
    _fitparticle(config().fitparticle()),
    _fdir(config().fitdirection()),
    _doSingleOutput(config().doSingleOutput()),
    _maxEDepAvg(config().maxEDepAvg()),
    _hfinder(config().hfinder()){
      consumes<ComboHitCollection>(_shLabel);
      consumes<TimeClusterCollection>(_timeclLabel);

      // Initialize the list of helicities to consider
      std::vector<int> helvals = config().Helicities();
      for(auto hv : helvals) {
          Helicity hel(hv);
          _hels.push_back(hel);
      }

      // In single output, put all helices by all helicities into one helix collection, otherwise produce a collection per helicity
      if (_doSingleOutput){
        produces<HelixSeedCollection>();
      } else {
        for(auto hel : _hels) {
          produces<HelixSeedCollection>(Helicity::name(hel));
        }
      }

      _tpart = (TrkParticle::type)_fitparticle;
//-----------------------------------------------------------------------------
// provide for interactive diagnostics
//-----------------------------------------------------------------------------

      if (_debugLevel != 0) _printfreq = 1;

      if (_diagLevel != 0) _hmanager = art::make_tool<ModuleHistToolBase>(config().diagPlugin,"diagPlugin");
      else                 _hmanager = std::make_unique<ModuleHistToolBase>();
    }

//-----------------------------------------------------------------------------
// destructor
//-----------------------------------------------------------------------------
  CalHelixFinder::~CalHelixFinder() {
  }

//-----------------------------------------------------------------------------
  void CalHelixFinder::beginJob(){
    if(_diagLevel != 0) { // only book histograms if running diagnostics
      art::ServiceHandle<art::TFileService> tfs;
      _hmanager->bookHistograms(tfs);
    }
  }

//-----------------------------------------------------------------------------
  void CalHelixFinder::beginRun(art::Run& ) {
    mu2e::GeomHandle<mu2e::BFieldManager> bfmgr;
    mu2e::GeomHandle<mu2e::DetectorSystem> det;
    Hep3Vector vpoint_mu2e = det->toMu2e(Hep3Vector(0.0,0.0,0.0));
    _bz0 = bfmgr->getBField(vpoint_mu2e).z();

    mu2e::GeomHandle<mu2e::Tracker> th;
    _tracker = th.get();

    mu2e::GeomHandle<mu2e::Calorimeter> ch;
    _calorimeter = ch.get();

    _hfinder.setTracker    (_tracker);
    _hfinder.setCalorimeter(_calorimeter);

    ChannelID cx, co;
    int       nPlanesPerStation(2);
    for (size_t ipl=0; ipl<_tracker->nPlanes(); ipl++) {
      const Plane*  pln = &_tracker->getPlane(ipl);
      for (size_t ipn=0; ipn<pln->nPanels(); ipn++) {
        const Panel* panel = &pln->getPanel(ipn);
        int face;
        if (panel->id().getPanel() % 2 == 0) face = 0;
        else                                 face = 1;
        cx.Station = ipl/nPlanesPerStation;//ist;
        cx.Plane   = ipl % nPlanesPerStation;
        cx.Face    = face;
        cx.Panel   = ipn;
        //            cx.Layer   = il;
        _hfResult.orderID (&cx, &co);
        int os = co.Station;
        int of = co.Face;
        int op = co.Panel;

        int       stationId = os;
        int       faceId    = of + stationId*StrawId::_nfaces*FaceZ_t::kNPlanesPerStation;
        _hfResult._zFace[faceId] = (panel->getStraw(0).getMidPoint().z()+panel->getStraw(1).getMidPoint().z())/2.;
        //-----------------------------------------------------------------------------
        // panel caches phi of its center and the z
        //-----------------------------------------------------------------------------
        _hfResult._phiPanel[faceId*FaceZ_t::kNPanels + op] = TVector2::Phi_0_2pi(polyAtan2(panel->straw0MidPoint().y(),panel->straw0MidPoint().x()));
      }
    }

    if (_debugLevel > 10){
      printf("//----------------------------------------------//\n");
      printf("//     Face      Panel       PHI       Z        //\n");
      printf("//----------------------------------------------//\n");

      for (int f=0; f<StrawId::_ntotalfaces; ++f){
        float z  =_hfResult._zFace[f];
        for (int p=0; p<FaceZ_t::kNPanels; ++p){
          float  phi = _hfResult._phiPanel[f*FaceZ_t::kNPanels + p];
          printf("//  %5i      %5i     %5.3f    %10.3f //\n", f, p, phi, z);
        }
      }
      printf("//----------------------------------//\n");

    }

  }

//-----------------------------------------------------------------------------
// find the input data objects
//-----------------------------------------------------------------------------
  bool CalHelixFinder::findData(const art::Event& evt) {

    if (evt.getByLabel(_shLabel, _strawhitsH)) {
      _chcol = _strawhitsH.product();
    }
    else {
      _chcol  = nullptr;
      printf(" >>> ERROR in CalHelixFinder::findData: StrawHitCollection with label=%s not found.\n",
             _shLabel.data());
    }

    if (evt.getByLabel(_timeclLabel, _timeclcolH)) {
      _timeclcol = _timeclcolH.product();
    }
    else {
      _timeclcol = nullptr;
      printf(" >>> ERROR in CalHelixFinder::findData: TimeClusterCollection with label=%s not found.\n",
             _timeclLabel.data());
    }
//-----------------------------------------------------------------------------
// done
//-----------------------------------------------------------------------------
   return (_chcol != nullptr) && (_timeclcol != nullptr);
  }

//-----------------------------------------------------------------------------
// event entry point
//-----------------------------------------------------------------------------
  void CalHelixFinder::produce(art::Event& event ) {
    const char*             oname = "CalHelixFinder::filter";

                                        // diagnostic info
    _data.event     = &event;
    _iev            = event.id().event();
    int   nGoodTClusterHits(0);

    if ((_debugLevel > 0) && (_iev%_printfreq) == 0) printf("[%s] : START event number %8i\n", oname,_iev);

    std::map<Helicity,std::unique_ptr<HelixSeedCollection>> helcols;
    int counter(0);
    if (!_doSingleOutput)  {
      for( auto const& hel : _hels) {
        helcols[hel] = std::unique_ptr<HelixSeedCollection>(new HelixSeedCollection());
        _data.nseeds [counter] = 0;
        ++counter;
      }
    } else {
      helcols[0] = std::unique_ptr<HelixSeedCollection>(new HelixSeedCollection());
      _data.nseeds [counter] = 0;
    }

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

    _data.nTimePeaks  = _timeclcol->size();
    for (int ipeak=0; ipeak<_data.nTimePeaks; ipeak++) {
      const TimeCluster* tc = &_timeclcol->at(ipeak);
      if (!tc->hasCaloCluster() && !_hfinder._procAllTCs)    continue;
      nGoodTClusterHits     = goodHitsTimeCluster(tc);
      if ( nGoodTClusterHits < _minNHitsTimeCluster)         continue;

      std::vector<HelixSeed>          helix_seed_vec;

//-----------------------------------------------------------------------------
// create track definitions for the helix fit from this initial information
// track fitting objects for this peak
//-----------------------------------------------------------------------------
      _hfResult.clearTempVariables();

      _hfResult._timeCluster    = tc;
      _hfResult._timeClusterPtr = art::Ptr<mu2e::TimeCluster>(_timeclcolH,ipeak);
//-----------------------------------------------------------------------------
// fill the face-order hits collector
//-----------------------------------------------------------------------------
      _hfinder.fillFaceOrderedHits(_hfResult);
//-----------------------------------------------------------------------------
// Step 1: now loop over the two possible helicities.
//         Find initial helical approximation of a track for both hypothesis
//-----------------------------------------------------------------------------
      std::vector<float> nHitsRatio_vec;
      for (size_t i=0; i<_hels.size(); ++i){
//-----------------------------------------------------------------------------
// create track definitions for the helix fit from this initial information
// track fitting objects for this peak
//-----------------------------------------------------------------------------
        CalHelixFinderData tmpResult(_hfResult);
        tmpResult.clearHelixInfo();

        tmpResult._helicity       = _hels[i];

        int rc = _hfinder.findHelix(tmpResult);

        if (!rc)                         continue;
        if(!tmpResult.helix()) {
          std::cout << "[CalHelixFinder::" << __func__ << "] " << event.id()
               << " Helix found but nullptr returned!!\n";
          continue;
        }
        HelixSeed     tmp_helix_seed;

        initHelixSeed(tmp_helix_seed, tmpResult);
        if (_diagLevel > 0) {
          nHitsRatio_vec.push_back(tmpResult._diag.nHitsRatio);
        }
        if (tmp_helix_seed._eDepAvg > _maxEDepAvg) continue;
        helix_seed_vec.push_back(tmp_helix_seed);
      }

      if (helix_seed_vec.size() == 0)                       continue;

//-----------------------------------------------------------------------------
// now select the best helix to avoid duplicates
//-----------------------------------------------------------------------------
      int    index_best(-1);
      pickBestHelix(helix_seed_vec, index_best);

//-----------------------------------------------------------------------------
// fill seed information
//-----------------------------------------------------------------------------
      if ( (index_best>=0) && (index_best < 2) ){
        Helicity              hel_best = helix_seed_vec[index_best]._helix._helicity;
        if (_doSingleOutput) {
          hel_best = 0;
        }
        HelixSeedCollection*  hcol     = helcols[hel_best].get();
        helix_seed_vec[index_best]._status.merge(TrkFitFlag::helixOK);
        hcol->push_back(helix_seed_vec[index_best]);
      } else if (index_best == 2){//both helices need to be saved

        for (unsigned k=0; k<_hels.size(); ++k){
          helix_seed_vec[k]._status.merge(TrkFitFlag::helixOK);
          Helicity              hel_best = helix_seed_vec[k]._helix._helicity;
          if (_doSingleOutput) {
            hel_best = 0;
          }
          HelixSeedCollection*  hcol     = helcols[hel_best].get();
          hcol->push_back(helix_seed_vec[k]);
        }
      }

      if (_diagLevel > 0) {
        fillDiagnosticInfo(helix_seed_vec, nHitsRatio_vec,
                           index_best, nGoodTClusterHits);
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
    int    nseeds(0);
    if (_doSingleOutput) {
      nseeds += helcols[0]->size();
      event.put(std::move(helcols[0]));
    }else    {
      for(auto const& hel : _hels ) {
        nseeds += helcols[hel]->size();
        event.put(std::move(helcols[hel]),Helicity::name(hel));
      }
    }
 }

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  void CalHelixFinder::endJob() {
  }

//--------------------------------------------------------------------------------
// set helix parameters
//-----------------------------------------------------------------------------
  void CalHelixFinder::initHelixSeed(HelixSeed& HelSeed, CalHelixFinderData& HfResult) {

    HelixTraj* hel = HfResult.helix();

    double   helixRadius   = 1./fabs(hel->omega());
    double   impactParam   = hel->d0();
    double   phi0          = hel->phi0();

    double   x0            = -(1/hel->omega()+impactParam)*sin(phi0);
    double   y0            =  (1/hel->omega()+impactParam)*cos(phi0);

    double   dfdz          = 1./hel->tanDip()/helixRadius;
    double   z0            = hel->z0();
                                        // center of the helix in the transverse plane
    Hep3Vector center(x0, y0, 0);
                                        //define the reconstructed helix parameters
    HelSeed._helix._rcent    = center.perp();
    HelSeed._helix._fcent    = center.phi();
    HelSeed._helix._radius   = helixRadius;
    HelSeed._helix._lambda   = 1./dfdz*_hfinder._dfdzsign;

    HelSeed._helix._fz0      = phi0 - M_PI/2.*_hfinder._dfdzsign -z0*hel->omega()/hel->tanDip() ;

    HelSeed._helix._helicity = HfResult._helicity;
    HelSeed._status.merge(TrkFitFlag::CPRHelix);

    //include also the values of the chi2d
    HelSeed._helix._chi2dXY   = HfResult._sxy.chi2DofCircle();
    HelSeed._helix._chi2dZPhi = HfResult._szphi.chi2DofLine();

                                        //now evaluate the helix T0 using the calorimeter cluster
    double   mm2MeV        = (3/10.)*_bz0;
    double   tandip        = hel->tanDip();
    double   mom           = helixRadius*mm2MeV/std::cos( std::atan(tandip));
    double   beta          = _tpart.beta(mom);

    double hel_t0;
    double pitchAngle = M_PI/2. - std::atan(tandip);

    HelSeed._timeCluster   = HfResult._timeClusterPtr;

    // cluster hits assigned to the reconsturcted Helix

    int nhits = HfResult.nGoodHits();

    HelSeed._hhits.setParent(_chcol->parent());

    for (int i=0; i<nhits; ++i){
      unsigned        hitId   = HfResult._goodhits[i];
      ComboHit*       hit     = &HfResult._chHitsToProcess[hitId];//panelz->_chHitsToProcess.at(hitInfo->panelHitIndex);
      ComboHit        hhit(*hit);
      HelSeed._hhits.push_back(hhit);
    }

    if(HfResult._timeClusterPtr->hasCaloCluster()) {
      CLHEP::Hep3Vector gpos,tpos;

      gpos   = _hfinder._calorimeter->geomUtil().diskToMu2e(HfResult._timeClusterPtr->caloCluster()->diskID(),
                                                                        HfResult._timeClusterPtr->caloCluster()->cog3Vector());
      tpos   = _hfinder._calorimeter->geomUtil().mu2eToTracker(gpos);

      hel_t0 = HfResult._timeClusterPtr->caloCluster()->time() - (tpos.z() - z0)/sin(pitchAngle)/(beta*CLHEP::c_light);

    }
    else {

      double tSum = 0.;

      size_t nHits = HelSeed.hits().size();

      if(nHits == 0) {
        hel_t0 = 0.;
      }
      else {

        for(size_t i = 0; i<nHits; i++) {

          const ComboHit& comboHit = HelSeed.hits()[i];
          double tCorrected = comboHit.correctedTime() - (comboHit.pos().z()-z0)/std::sin(pitchAngle)/(beta*CLHEP::c_light);
          tSum += tCorrected;

        }
        hel_t0 = tSum/((double) nHits);
      }

    }

    HelSeed._t0            = TrkT0(hel_t0, 0.1); //dummy error on T0 FIXME!

    HelSeed._eDepAvg = HelSeed._hhits.eDepAvg();


    //now set the HelixRecoDir
    HelixTool ht(&HelSeed, _tracker);
    float     slope(0), slopeErr(0), chi2ndof(0);
    ht.dirOfProp(slope, slopeErr, chi2ndof);
    HelSeed._recoDir._slope    = slope;
    HelSeed._recoDir._slopeErr = slopeErr;
    HelSeed._recoDir._chi2ndof = chi2ndof;
  }

//-----------------------------------------------------------------------------
  int CalHelixFinder::initHelixFinderData(CalHelixFinderData&                Data,
                                          const TrkParticle&                 TPart,
                                          const TrkFitDirection&             FDir,
                                          const ComboHitCollection*          ComboCollection ) {
    Data._fit         = TrkErrCode::fail;
    Data._tpart       = TPart;
    Data._fdir        = FDir;

    Data._chcol       = ComboCollection;

    Data._radius      = -1.0;
    Data._dfdz        = 0.;
    Data._fz0         = 0.;

    return 0;
  }

//-----------------------------------------------------------------------------
  int  CalHelixFinder::goodHitsTimeCluster(const TimeCluster* TCluster){
    int   nhits         = TCluster->nhits();
    int   ngoodhits(0);
    for (int i=0; i<nhits; ++i){
      int          index   = TCluster->hits().at(i);
      ComboHit     sh      = _chcol ->at(index);
      StrawHitFlag flag    = sh.flag();
      int          bkg_hit = flag.hasAnyProperty(StrawHitFlag::bkg);
      if (bkg_hit)                              continue;

      ngoodhits += sh.nStrawHits();
    }

    return ngoodhits;
  }

//--------------------------------------------------------------------------------
// function to select the best Helix among the results of the two helicity hypo
//--------------------------------------------------------------------------------
  void  CalHelixFinder::pickBestHelix(std::vector<HelixSeed>& HelVec, int &Index_best){
    if (HelVec.size() == 1) {
      Index_best = 0;
      return;
    }

//-----------------------------------------------------------------------------
// at Mu2e, 2 helices with different helicity could be duplicates of each other
//-----------------------------------------------------------------------------

    const HelixSeed           *h1, *h2;
    const ComboHitCollection  *tlist, *clist;
    int                        nh1, nh2, natc(0);

    constexpr float overlap_threshold(0.5f); // 50% overlap to match helices

    // Helix 1
    h1     = &HelVec[0];
    tlist  = &h1->hits();
    nh1    = tlist->size();

    // Helix 2
    h2     = &HelVec[1];
    clist  = &h2->hits();
    nh2    = clist->size();

//-----------------------------------------------------------------------------
// check the number of common hits
//-----------------------------------------------------------------------------
    std::unordered_set<StrawHitIndex> hits_1; // unordered set for O(1) lookups
    hits_1.reserve(nh1);

    // Fill the set
    for(int k = 0; k < nh1; ++k) {
      hits_1.insert(tlist->at(k).index());
    }

    // Count overlapping hits
    for(int l = 0; l < nh2; l++) {
      if(hits_1.contains(clist->at(l).index())) ++natc;
    }


    // Check if at least one shares half of its hits with the other helix
    if ((natc > nh1*overlap_threshold) || (natc > nh2*overlap_threshold)) {

 //-----------------------------------------------------------------------------
 // pick the helix with the largest number of hits
 //-----------------------------------------------------------------------------
      if      (nh2 > nh1) Index_best = 1; // helix two has more hits
      else if (nh1 > nh2) Index_best = 0; // helix one has more hits
      else                Index_best = (h1->helix().chi2dZPhi() < h2->helix().chi2dZPhi()) ? 0 : 1; // equal hits --> pick the best chi^2
    } else {
//-----------------------------------------------------------------------------
// this is the case where we consider the two helices independent, so we want
// to store both
//-----------------------------------------------------------------------------
      Index_best  = 2;
    }


  }


//-----------------------------------------------------------------------------
// fill diagnostic information
//-----------------------------------------------------------------------------
  void CalHelixFinder::fillDiagnosticInfo(const std::vector<HelixSeed>& helix_seed_vec,
                                          const std::vector<float>& nHitsRatio_vec,
                                          int index_best, const int nGoodTClusterHits) {
    const char*     oname = "CalHelixFinder::fillDiagnosticsInfo";
    int             nhitsMin(15);
    double          mm2MeV = (3/10.)*_bz0;

    int loc = _data.nseeds[0];
    if (loc < _data.maxSeeds()) {
      if (index_best == 2){
        index_best = 0;
      }
      int nhits          = helix_seed_vec[index_best]._hhits.size();
      _data.ntclhits[loc]= nGoodTClusterHits;
      _data.nhits[loc]   = nhits;
      _data.radius[loc]  = helix_seed_vec[index_best].helix().radius();
      _data.pT[loc]      = mm2MeV*_data.radius[loc];
      _data.p[loc]       = _data.pT[loc]/std::cos( std::atan(helix_seed_vec[index_best].helix().lambda()/_data.radius[loc]));

      _data.chi2XY[loc]   = _hfResult._sxy.chi2DofCircle();
      _data.chi2ZPhi[loc] = _hfResult._szphi.chi2DofLine();

      _data.nHitsRatio[loc] = nHitsRatio_vec[index_best];
      _data.eDepAvg[loc] = helix_seed_vec[index_best]._eDepAvg;

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
      for (unsigned i=0; i<_hfResult._hitsUsed.size(); ++i){
        if (_hfResult._hitsUsed[i] != 1)           continue;
        ++counter;
      }
    }
    else {
      printf(" N(seeds) > %i, IGNORE SEED\n",_data.maxSeeds());
    }
  }
}

using mu2e::CalHelixFinder;
DEFINE_ART_MODULE(CalHelixFinder)
