///////////////////////////////////////////////////////////////////////////////
// Calorimeter-driven track finding
// Search for clusters of strahits using the calorimeter cluster
// It passes TimeClusters to CalPatRecNew
// P.Murat, G.Pezzullo
// try to order routines alphabetically
///////////////////////////////////////////////////////////////////////////////
#include "fhiclcpp/ParameterSet.h"

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

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/algorithm/string.hpp>

#include "TVector2.h"

using namespace std;
using namespace boost::accumulators;
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
    _ccmLabel        (pset.get<string>         ("caloClusterModuleLabel"         )),
    _hsel            (pset.get<vector<string> >("HelixFitSelectionBits"          )),
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

    if (_diagLevel  != 0) _hmanager = art::make_tool<CprModuleHistBase>(pset.get<fhicl::ParameterSet>("histograms"));
    else                  _hmanager = std::make_unique<CprModuleHistBase>();
  }

//-----------------------------------------------------------------------------
// destructor
//-----------------------------------------------------------------------------
  CalTimePeakFinder::~CalTimePeakFinder() {}

//-----------------------------------------------------------------------------
  void CalTimePeakFinder::beginJob(){
    art::ServiceHandle<art::TFileService> tfs;
    _hmanager->bookHistograms(tfs,&_hist);
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
// find the input data objects
//-----------------------------------------------------------------------------
  bool CalTimePeakFinder::findData(const art::Event& evt) {

    //    art::Handle<mu2e::StrawHitCollection> strawhitsH;
    if (evt.getByLabel(_shLabel, _strawhitsH)) {
      _shcol = _strawhitsH.product();
    }
    else {
      _shcol  = 0;
      printf(" >>> ERROR in CalTimePeakFinder::findData: StrawHitCollection with label=%s not found.\n",
             _shLabel.data());
    }

    art::Handle<mu2e::StrawHitFlagCollection> shflagH;
    if (evt.getByLabel(_shfLabel,shflagH)) {
      _shfcol = shflagH.product();
    }
    else {
      _shfcol = 0;
      printf(" >>> ERROR in CalTimePeakFinder::findData: StrawHitFlagCollection with label=%s not found.\n",
             _shfLabel.data());
    }

    if (evt.getByLabel(_ccmLabel, _ccH)) {
      _ccCollection = _ccH.product();
    }
    else {
      _ccCollection = 0;
      printf(" >>> ERROR in CalTimePeakFinder::findData: CaloClusterCollection with label=%s not found.\n",
             _ccmLabel.data());
    }

 
//-----------------------------------------------------------------------------
// done
//-----------------------------------------------------------------------------
    return (_shcol != 0) && (_shfcol != 0) && (_ccCollection != 0);
  }

//-----------------------------------------------------------------------------
// event entry point
//-----------------------------------------------------------------------------
//  void CalTimePeakFinder::produce(art::Event& event ) {
  bool CalTimePeakFinder::filter(art::Event& event ) {
    const char*               oname = "CalTimePeakFinder::filter";

                                        // event printout
    _iev     = event.id().event();

    if ((_iev%_printfreq) == 0) printf("[%s] : START event number %8i\n", oname,_iev);

    _tpeaks = new CalTimePeakCollection;
    
    unique_ptr<TimeClusterCollection>  outseeds(new TimeClusterCollection);
    unique_ptr<CalTimePeakCollection>  tpeaks  (_tpeaks);
//-----------------------------------------------------------------------------
// find the time peaks in the time spectrum of selected hits.
//-----------------------------------------------------------------------------
    bool ok = findData(event);

    if (ok) findTimePeaks(_tpeaks, *outseeds);
    else {
      printf("%s ERROR: No straw hits found\n",oname);
    }
//--------------------------------------------------------------------------------    
// fill diagnostic if needed
//--------------------------------------------------------------------------------
    if (_diagLevel > 0) {

      int   nseeds = outseeds->size();
      
      _data.nseeds[0] = nseeds;
      _data.nseeds[1] = 0;
      _data.minNHits  = _minNHits;
      
      for (int i=0; i<nseeds; ++i) {
	TimeCluster* tmpseed = &outseeds->at(i);
	_data.cl             = tmpseed->caloCluster().get();
	_data.timeCluster    = tmpseed;
	int nhits            = tmpseed->hits().size();

	if (nhits >= _minNHits) _data.nseeds[1] += 1;
//-----------------------------------------------------------------------------
// fill timepeak-level histograms
//-----------------------------------------------------------------------------
	_hmanager->fillHistograms(1,&_data,&_hist);
      }
//-----------------------------------------------------------------------------
// fill event-level histograms : so far, for nseeds
//-----------------------------------------------------------------------------
      _hmanager->fillHistograms(0,&_data,&_hist);
    }
//-----------------------------------------------------------------------------
// put reconstructed tracks into the event record
//-----------------------------------------------------------------------------
    int   nseeds = outseeds->size();
    
    event.put(std::move(outseeds));
    event.put(std::move(tpeaks  ));
    
    if (_useAsFilter == 1) return (nseeds > 0) ? true : false ; 
    else                   return true;

  }

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  void CalTimePeakFinder::endJob(){ }

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  void CalTimePeakFinder::findTimePeaks(CalTimePeakCollection* TimePeakColl, TimeClusterCollection& OutSeeds) {

    const char* oname = "CalTimePeakFinder::findTimePeaks";

    int                 ncl, nsh;
    double              time, dt, tof, zstraw, cl_time, stime;
    double              xcl, ycl, zcl/*, dz_cl*/;
    const CaloCluster*  cl;
    const StrawHit*     hit;
    const Straw*        straw;
    Hep3Vector          gpos, tpos;

    StrawHitFlag        energyFlag(StrawHitFlag::energysel);
    StrawHitFlag        timeFlag  (StrawHitFlag::timesel);
    StrawHitFlag        radiusFlag(StrawHitFlag::radsel);
    StrawHitFlag        deltaRayFlag(StrawHitFlag::delta);
    StrawHitFlag        isolatedFlag(StrawHitFlag::isolated);
//-----------------------------------------------------------------------------
// Loop over calorimeter clusters
//-----------------------------------------------------------------------------
    nsh   = _shcol->size();
    ncl   = _ccCollection->size();

    for (int ic=0; ic<ncl; ic++) {
      cl      = &_ccCollection->at(ic);

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

          //    dz_cl   = zcl; // -_tracker->z0();
          // create time peak
          CalTimePeak tpeak(cl,xcl,ycl,zcl);

          tpeak._shcol  = _shcol;
          tpeak._shfcol = _shfcol;
          tpeak._tmin   = cl_time+_mindt;
          tpeak._tmax   = cl_time+_maxdt;
//-----------------------------------------------------------------------------
// record hits in time with each peak, and accept them if they have a minimum # of hits
//-----------------------------------------------------------------------------
          stime = 0;
          mu2e::StrawHitFlag flag;
          // int   nhitsTimeWindow(0), nhitsHasTime(0), nhitsHasEnergy(0), nhitsHasRadius(0),
          //       nhitsNoDelta(0), nhitsNoIsolated(0);

          double meanDriftTime = 1.25/0.06;// half straw tube radius / drift velocity
	  //          int    gen_index, sim_id, vol_id;

          for(int istr=0; istr<nsh;++istr) {
            flag = _shfcol->at(istr);

            int hit_has_all_properties = flag.hasAllProperties(_hsel);
            int bgr_hit                = flag.hasAnyProperty(_bkgsel);

            // int hit_has_energy         = flag.hasAllProperties(energyFlag);
            // int hit_has_time           = flag.hasAllProperties(timeFlag);
            // int hit_has_radius         = flag.hasAllProperties(radiusFlag);

            // int deltaRay_hit           = flag.hasAnyProperty(deltaRayFlag);
            // int isolated_hit           = flag.hasAnyProperty(isolatedFlag);

            hit    = &_shcol->at(istr);
            time   = hit->time();
            straw  = &_tracker->getStraw(hit->strawIndex());
            zstraw = straw->getMidPoint().z();
//-----------------------------------------------------------------------------
// estimate time-of-flight and calculate residual between the predicted and the hit times
//-----------------------------------------------------------------------------
            tof = (zcl-zstraw)/sin(_pitchAngle)/CLHEP::c_light;
            dt  = cl_time-(time+tof-meanDriftTime);
//-----------------------------------------------------------------------------
// fill some diag histograms
//-----------------------------------------------------------------------------
            if ((dt < _maxdt) && (dt >= _mindt)) {

              if (hit_has_all_properties && !bgr_hit) {
                tpeak._index.push_back(istr);
                stime += time;
              }
	      else if (_debugLevel > 0) {
//-----------------------------------------------------------------------------
// print diagnostics on rejected hits
//-----------------------------------------------------------------------------
		printf("[%s] rejected hit: index: %5i flag: %10s  time:  %8.3f   dt: %8.3f energy: %8.5f\n",
		       oname, istr, flag.hex().data(), hit->time(), hit->dt(), hit->energyDep());
	      }
            }
          }

          tpeak._tpeak = stime/(tpeak.NHits()+1.e-12);

          if (tpeak.NHits() > _minNHits) {
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
