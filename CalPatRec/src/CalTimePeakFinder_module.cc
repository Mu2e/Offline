///////////////////////////////////////////////////////////////////////////////
// Calorimeter-driven track finding
// Search for clusters of strahits using the calorimeter cluster
// It passes TimeClusters to CalPatRecNew
// P.Murat, G.Pezzullo
// try to order routines alphabetically
// *FIXME* : need to use the assumed particle velocity instead of the speed of light
///////////////////////////////////////////////////////////////////////////////
#include "fhiclcpp/ParameterSet.h"

#include "CalPatRec/inc/ModuleHistToolBase.hh"
#include "CalPatRec/inc/CalTimePeakFinder_module.hh"

// framework
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Utilities/make_tool.h"

// conditions
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"

#include "TTrackerGeom/inc/TTracker.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"

#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"

#include "RecoDataProducts/inc/StrawHitIndex.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"

using namespace std;

using CLHEP::HepVector;
using CLHEP::HepSymMatrix;
using CLHEP::Hep3Vector;

namespace mu2e {
//-----------------------------------------------------------------------------
// comparison functor for sorting by Z(wire)
//-----------------------------------------------------------------------------
  struct straw_zcomp : public binary_function<StrawHitIndex,StrawHitIndex,bool> {
    bool operator()(StrawHitIndex const& h1, StrawHitIndex const& h2) {

      mu2e::GeomHandle<mu2e::TTracker> handle;
      const TTracker* t = handle.get();
      const Straw* s1 = &t->getStraw(StrawIndex(h1));
      const Straw* s2 = &t->getStraw(StrawIndex(h2));

      return s1->getMidPoint().z() < s2->getMidPoint().z();
    }
  }; // a semicolumn here is required

//-----------------------------------------------------------------------------
// module constructor, parameter defaults are defiend in CalPatRec/fcl/prolog.fcl
//-----------------------------------------------------------------------------
  CalTimePeakFinder::CalTimePeakFinder(fhicl::ParameterSet const& pset) :
    _diagLevel       (pset.get<int>            ("diagLevel"                      )),
    _debugLevel      (pset.get<int>            ("debugLevel"                     )),
    _printfreq       (pset.get<int>            ("printFrequency"                 )),
    _useAsFilter     (pset.get<int>            ("useAsFilter"                    )),    
    _shLabel         (pset.get<string>         ("StrawHitCollectionLabel"        )),
    _shfLabel        (pset.get<string>         ("StrawHitFlagCollectionLabel"    )),
    _shpLabel        (pset.get<string>         ("StrawHitPositionCollectionLabel")),
    _ccmLabel        (pset.get<string>         ("caloClusterModuleLabel"         )),
    _hsel            (pset.get<vector<string> >("HitSelectionBits"               )),
    _bkgsel          (pset.get<vector<string> >("BackgroundSelectionBits"        )),
    _mindt           (pset.get<double>         ("DtMin"                          )),
    _maxdt           (pset.get<double>         ("DtMax"                          )),
    _minNHits        (pset.get<int>            ("MinNHits"                       )),
    _minClusterEnergy(pset.get<double>         ("minClusterEnergy"               )),
    _minClusterSize  (pset.get<int>            ("minClusterSize"                 )),
    _minClusterTime  (pset.get<double>         ("minClusterTime"                 )),
    _pitchAngle      (pset.get<double>         ("pitchAngle"                     )),
    _dtoffset        (pset.get<double>         ("dtOffset"                       ))
  {
    produces<TimeClusterCollection>();
    produces<CalTimePeakCollection>();

    if (_debugLevel != 0) _printfreq = 1;

    if (_diagLevel  != 0) _hmanager = art::make_tool<ModuleHistToolBase>(pset.get<fhicl::ParameterSet>("diagPlugin"));
    else                  _hmanager = std::make_unique<ModuleHistToolBase>();

    _data.minClusterEnergy =  _minClusterEnergy;
    _data.minNHits         =  _minNHits;
  }

//-----------------------------------------------------------------------------
// destructor
//-----------------------------------------------------------------------------
  CalTimePeakFinder::~CalTimePeakFinder() {}

//-----------------------------------------------------------------------------
  void CalTimePeakFinder::beginJob(){
    if (_diagLevel > 0) {
      art::ServiceHandle<art::TFileService> tfs;
      _hmanager->bookHistograms(tfs);
    }
  }

//-----------------------------------------------------------------------------
  bool CalTimePeakFinder::beginRun(art::Run& ) {
    mu2e::GeomHandle<mu2e::TTracker> th;
    _tracker = th.get();

    mu2e::GeomHandle<mu2e::Calorimeter> ch;
    _calorimeter = ch.get();
					// calibrations

    mu2e::ConditionsHandle<TrackerCalibrations> tcal("ignored");
    _trackerCalib = tcal.operator ->();
    
    return true;
  }

//-----------------------------------------------------------------------------
// find input things
//-----------------------------------------------------------------------------
  bool CalTimePeakFinder::findData(const art::Event& evt) {

    //    art::Handle<mu2e::StrawHitCollection> shcolH;
    auto shcolH = evt.getValidHandle<mu2e::StrawHitCollection>(_shLabel);
    if (shcolH.product() != 0){
      //    if (evt.getByLabel(_shLabel, shcolH)) {
      _data.shcol = shcolH.product();
    }
    else {
      _data.shcol  = 0;
      printf(" >>> ERROR in CalTimePeakFinder::findData: StrawHitCollection with label=%s not found.\n",
             _shLabel.data());
    }


    //    art::Handle<mu2e::StrawHitPositionCollection> shposH;
    auto shposH = evt.getValidHandle<mu2e::StrawHitPositionCollection>(_shpLabel);
    if (shcolH.product() != 0){
      //    if (evt.getByLabel(_shpLabel,shposH)) {
      _data.shpcol = shposH.product();
    }
    else {
      _data.shpcol = 0;
      printf(" >>> ERROR in CalHelixFinder::findData: StrawHitPositionCollection with label=%s not found.\n",
             _shpLabel.data());
    }

    //    art::Handle<mu2e::StrawHitFlagCollection> shflagH;
    // auto shflagH = evt.getValidHandle<mu2e::StrawHitFlagCollection>(_shfLabel);
    // if (shflagH.product() != 0){
    // //    if (evt.getByLabel(_shfLabel,shflagH)) {
    //   _data.shfcol = shflagH.product();
    // }
    // else {
    //   _data.shfcol = 0;
    //   printf(" >>> ERROR in CalTimePeakFinder::findData: StrawHitFlagCollection with label=%s not found.\n",
    //          _shfLabel.data());
    // }

    if (evt.getByLabel(_ccmLabel, _ccH)) {
      _data.ccCollection = _ccH.product();
    }
    else {
      _data.ccCollection = 0;
      printf(" >>> ERROR in CalTimePeakFinder::findData: CaloClusterCollection with label=%s not found.\n",
             _ccmLabel.data());
    }

    return (_data.shcol != 0) && /*(_data.shfcol != 0) && */(_data.ccCollection != 0);
  }

//-----------------------------------------------------------------------------
// event entry point
//-----------------------------------------------------------------------------
  bool CalTimePeakFinder::filter(art::Event& event) {
    const char*               oname = "CalTimePeakFinder::filter";

                                        // event printout
    _iev     = event.id().event();
    if ((_iev%_printfreq) == 0) printf("[%s] : START event number %8i\n", oname,_iev);

    _data._event = &event;
    _data._tpeaks = new CalTimePeakCollection;

    unique_ptr<CalTimePeakCollection>  tpeaks  (_data._tpeaks);
    unique_ptr<TimeClusterCollection>  outseeds(new TimeClusterCollection);
    
    _data._outseeds = outseeds.get();

    bool ok = findData(event);

    if (ok) findTimePeaks(_data._tpeaks, *_data._outseeds);
    else    printf("%s ERROR: No straw hits found in event %i\n",oname,_iev);

    // diagnostics, if requested
    if (_diagLevel > 0) _hmanager->fillHistograms(&_data);
//-----------------------------------------------------------------------------
// put reconstructed tracks into the event record
//-----------------------------------------------------------------------------
    event.put(std::move(outseeds));
    event.put(std::move(tpeaks  ));
//-----------------------------------------------------------------------------
// filtering, if requested
//-----------------------------------------------------------------------------
    if (_useAsFilter == 0) return true;

    int nseeds = outseeds->size();
    return (nseeds > 0) ; 
  }
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  void CalTimePeakFinder::endJob(){ }

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  void CalTimePeakFinder::findTimePeaks(CalTimePeakCollection* TimePeakColl, TimeClusterCollection& OutSeeds) {

    //    const char* oname = "CalTimePeakFinder::findTimePeaks";

    int                 ncl, nsh;
    double              time, dt, tof, zstraw, cl_time, stime;
    double              xcl, ycl, zcl/*, dz_cl*/;
    const CaloCluster*  cl;
    const StrawHit*     hit;
    const Straw*        straw;
    Hep3Vector          gpos, tpos;

    //    using namespace boost::accumulators;
    static const double pi(M_PI);
    static const double twopi(2*pi);

    double              mphi(-9999.);

    StrawHitFlag        energyFlag(StrawHitFlag::energysel);
    StrawHitFlag        timeFlag  (StrawHitFlag::timesel);
    StrawHitFlag        radiusFlag(StrawHitFlag::radsel);
    StrawHitFlag        deltaRayFlag(StrawHitFlag::bkg);
    StrawHitFlag        isolatedFlag(StrawHitFlag::isolated);
//-----------------------------------------------------------------------------
// Loop over calorimeter clusters
//-----------------------------------------------------------------------------
    nsh   = _data.shcol->size();
    ncl   = _data.ccCollection->size();

    for (int ic=0; ic<ncl; ic++) {
      cl      = &_data.ccCollection->at(ic);

      if ( cl->energyDep() > _minClusterEnergy) {

        if ( (int(cl->size()) >= _minClusterSize) ) {

          cl_time = cl->time() + _dtoffset;
//-----------------------------------------------------------------------------
// convert cluster coordinates defined in the disk frame to the detector
// coordinate system
//-----------------------------------------------------------------------------
          gpos = _calorimeter->geomUtil().diskToMu2e(cl->diskId(),cl->cog3Vector());
          tpos = _calorimeter->geomUtil().mu2eToTracker(gpos);

          xcl     = tpos.x();
          ycl     = tpos.y();
          zcl     = tpos.z();
	  mphi    = atan2(ycl, xcl);

          //    dz_cl   = zcl; // -_tracker->z0();
          // create time peak
          CalTimePeak tpeak(cl,xcl,ycl,zcl);

          tpeak._shcol  = _data.shcol;
          tpeak._shfcol = NULL; // _data.shfcol;
          tpeak._tmin   = cl_time+_mindt;
          tpeak._tmax   = cl_time+_maxdt;
//-----------------------------------------------------------------------------
// record hits in time with each peak, and accept them if they have a minimum # of hits
//-----------------------------------------------------------------------------
          stime = 0;
          mu2e::StrawHitFlag flag;

          double meanDriftTime = 1.25/0.06;// half straw tube radius / drift velocity

          for(int istr=0; istr<nsh;++istr) {
	    //            flag = _data.shfcol->at(istr);

            // int hit_has_all_properties = flag.hasAllProperties(_hsel);
            // int bgr_hit                = flag.hasAnyProperty(_bkgsel);

            hit    = &_data.shcol->at(istr);
            time   = hit->time();
            straw  = &_tracker->getStraw(hit->strawIndex());
            zstraw = straw->getMidPoint().z();
//-----------------------------------------------------------------------------
// estimate time-of-flight and calculate residual between the predicted and the hit times
// 2017-03-31 P.M.: this assumes electron (e^- or e^+), not muon
//-----------------------------------------------------------------------------
            tof = (zcl-zstraw)/sin(_pitchAngle)/CLHEP::c_light;
            dt  = cl_time-(time+tof-meanDriftTime);
//--------------------------------------------------------------------------------
// check the angular distance from the calorimeter cluster
// 2017-11-17 Gianipez: this selection was present on CalHelixFinderAlg::filterDist
//--------------------------------------------------------------------------------
	    double dphi = _data.shpcol->at(istr).pos().phi() - mphi;
	    
	    if (dphi >  pi) dphi -= twopi;
	    if (dphi < -pi) dphi += twopi;

//-----------------------------------------------------------------------------
// fill some diag histograms
//-----------------------------------------------------------------------------
            if ((dt < _maxdt) && (dt >= _mindt) && (fabs(dphi) <= pi/2.) ) {

	      //              if (hit_has_all_properties && !bgr_hit) {
	      tpeak._index.push_back(istr);
	      stime += time;
		//              }
	      //else if (_debugLevel > 0) {
//-----------------------------------------------------------------------------
// print diagnostics on rejected hits
//-----------------------------------------------------------------------------
		// printf("[%s] rejected hit: index: %5i flag: %10s  time:  %8.3f   dt: %8.3f energy: %8.5f\n",
		//        oname, istr, flag.hex().data(), hit->time(), hit->dt(), hit->energyDep());
	      // }
            }
          }

          tpeak._tpeak = stime/(tpeak.NHits()+1.e-12);

          if (tpeak.NHits() >= _data.minNHits) {
	    TimePeakColl->push_back(tpeak);

					//fill seed information
	    TimeCluster    tmpseed;
	    initTimeCluster(tmpseed, tpeak, ic);
	    OutSeeds.push_back(tmpseed);
          }
        }
      }
    }
  }

//--------------------------------------------------------------------------------
  void CalTimePeakFinder::initTimeCluster(TimeCluster  &TrkSeed     , 
					  CalTimePeak  &TPeak       , 
					  int          &ClusterIndex) {
    
    int             shIndices = TPeak.NHits();

    for (int i=0; i<shIndices; ++i){
      size_t   hIndex = TPeak._index.at(i);
      TrkSeed._strawHitIdxs.push_back( StrawHitIndex( hIndex) );
    }
    
    const mu2e::CaloCluster *cluster = TPeak.Cluster();
    
                                //do we need to propagate the cluster time at z=0 at this stage?
    TrkSeed._t0               = TrkT0(cluster->time(), 0.1); //dummy value for errT0
    
    int               idisk   = cluster->diskId();
    CLHEP::Hep3Vector cp_mu2e = _calorimeter->geomUtil().diskToMu2e(idisk, cluster->cog3Vector());
    CLHEP::Hep3Vector cp_st   = _calorimeter->geomUtil().mu2eToTracker(cp_mu2e);
    
    
    TrkSeed._pos              = cp_st;
    TrkSeed._caloCluster      = art::Ptr<mu2e::CaloCluster>(_ccH, ClusterIndex);
  }


}

using mu2e::CalTimePeakFinder;
DEFINE_ART_MODULE(CalTimePeakFinder);
