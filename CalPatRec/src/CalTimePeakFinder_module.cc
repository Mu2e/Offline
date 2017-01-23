///////////////////////////////////////////////////////////////////////////////
// Calorimeter-driven track finding
// Search for clusters of strahits using the calorimeter cluster
// It passes TimeClusters to CalPatRecNew
// P.Murat, G.Pezzullo
// try to order routines alphabetically
///////////////////////////////////////////////////////////////////////////////
#include "fhiclcpp/ParameterSet.h"

#include "CalPatRec/inc/CalTimePeakFinder_module.hh"
#include "CalPatRec/inc/Ref.hh"

// framework
#include "art/Framework/Principal/Handle.h"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"

// conditions
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "CalPatRec/inc/KalFitResult.hh"
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
    _minnhits        (pset.get<int>            ("MinNHits"                       )),
    _minClusterEnergy(pset.get<double>         ("minClusterEnergy"               )),
    _minClusterSize  (pset.get<int>            ("minClusterSize"                 )),
    _minClusterTime  (pset.get<double>         ("minClusterTime"                 )),
    _pitchAngle      (pset.get<double>         ("pitchAngle"                     )),
    _dtoffset        (pset.get<double>         ("dtOffset"                       ))
  {
    produces<TimeClusterCollection>();
    produces<CalTimePeakCollection>();

  }

//-----------------------------------------------------------------------------
// destructor
//-----------------------------------------------------------------------------
  CalTimePeakFinder::~CalTimePeakFinder() {}

//-----------------------------------------------------------------------------
  void CalTimePeakFinder::beginJob(){

    if(_diagLevel > 0) bookHistograms();

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
//
//-----------------------------------------------------------------------------
  void CalTimePeakFinder::bookHistograms() {
    art::ServiceHandle<art::TFileService> tfs;

    art::TFileDirectory hf_dir = tfs->mkdir("TimePeak");
    
    _hist.nseeds[0]          = tfs->make<TH1F>("nseeds0"  , "number of track candidates: all events", 21, -0.5, 20.5);
    _hist.nseeds[1]          = tfs->make<TH1F>("nseeds1"  , "number of track candidates: nhits > 15", 21, -0.5, 20.5);
    _hist.timePeak.nhits     = hf_dir.make<TH1F>("nhits"  , "number of hits within a track candidate; nHits", 101, -0.5, 100.5);
    _hist.timePeak.energy[0] = hf_dir.make<TH1F>("energy0", "cluster energy; E [MeV]"                   , 400, 0., 200.);
    _hist.timePeak.energy[1] = hf_dir.make<TH1F>("energy1", "cluster energy, nhits > 15; E [MeV]"       , 400, 0., 200.);
    _hist.timePeak.time  [0] = hf_dir.make<TH1F>("time0"  , "cluster time; t [ns]"                      , 2800, 300., 1700);
    _hist.timePeak.time  [1] = hf_dir.make<TH1F>("time1"  , "cluster time, nhits > 15; t [ns]"          , 2800, 300., 1700);

    _hist.timePeak.nhitsvstime   = hf_dir.make<TH2F>("nhitsvstime","nhits vs time; N [#]; t [ns]"       , 100, 0, 100, 2800, 300., 1700);
    _hist.timePeak.nhitsvsenergy = hf_dir.make<TH2F>("nhitsvsenergy" ,"nhits vs energy; N [#]; E [MeV]" , 100, 0, 100, 400, 0, 200);

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
    const char*               oname = "CalTimePeakFinder::produce";

                                        // event printout
    _iev     = event.id().event();

    if ((_iev%_printfreq) == 0) printf("[%s] : START event number %8i\n", oname,_iev);

    _tpeaks = new CalTimePeakCollection;
    
    unique_ptr<TimeClusterCollection>  outseeds(new TimeClusterCollection);
    unique_ptr<CalTimePeakCollection>  tpeaks  (_tpeaks);

    
                                        // find the data
    if (!findData(event)) {
      printf("%s ERROR: No straw hits found, RETURN\n",oname);
                                                            goto END;
    }

//-----------------------------------------------------------------------------
// find the time peaks in the time spectrum of selected hits.
//-----------------------------------------------------------------------------
    findTimePeaks(_tpeaks, *outseeds);

//--------------------------------------------------------------------------------    
// fill diagnostic if needed
//--------------------------------------------------------------------------------
    if (_diagLevel > 0) {
      int   nseeds = outseeds->size();
      
      _hist.nseeds[0]->Fill(nseeds);
      
      double                     clTime(0), clEnergy(0), nseedsCut0(0);
      int                        nhits(0);
      TimeCluster               *tmpseed;
      const      CaloCluster    *cluster;
      
      for (int i=0; i<nseeds; ++i){
	tmpseed    = &outseeds->at(i);
	cluster    = tmpseed->caloCluster().get();

	clTime     = cluster->time();
	clEnergy   = cluster->energyDep();
	nhits      = tmpseed->hits().size();

	if (nhits >= 15) {
	  ++nseedsCut0;
	  _hist.timePeak.energy[1] ->Fill(clEnergy);
	  _hist.timePeak.time  [1] ->Fill(clTime);
	}
  	
	_hist.timePeak.energy[0] ->Fill(clEnergy);
	_hist.timePeak.time  [0] ->Fill(clTime);
	_hist.timePeak.nhits     ->Fill(nhits);

	_hist.timePeak.nhitsvstime   ->Fill(nhits, clEnergy);
	_hist.timePeak.nhitsvsenergy ->Fill(nhits, clTime);
	
      }

      _hist.nseeds[1]->Fill(nseedsCut0);

      
    }

//-----------------------------------------------------------------------------
// put reconstructed tracks into the event record
//-----------------------------------------------------------------------------
  END:;
    int   nseeds = outseeds->size();
    
    event.put(std::move(outseeds));
    event.put(std::move(tpeaks));
    
    if (_useAsFilter == 1) {
      if (nseeds > 0) {
	return true;
      } else{
	return false;
      }
    } else {
      return true;
    }


  }

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  void CalTimePeakFinder::endJob(){ }

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  void CalTimePeakFinder::findTimePeaks(CalTimePeakCollection* TimePeakColl, TimeClusterCollection& OutSeeds) {

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
          gpos = _calorimeter->geomUtil().mu2eToDiskFF(cl->diskId(),cl->cog3Vector());
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
            }
          }

          tpeak._tpeak = stime/(tpeak.NHits()+1.e-12);

          if (tpeak.NHits() > _minnhits)       {
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
  void CalTimePeakFinder::initTimeCluster(TimeCluster                           &TrkSeed     , 
					  CalTimePeak                           &TPeak       , 
					  int                                   &ClusterIndex){
    
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
