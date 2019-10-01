///////////////////////////////////////////////////////////////////////////////
// Calorimeter-driven track finding
// Search for clusters of strahits using the calorimeter cluster
// It passes TimeClusters to CalPatRecNew
// P.Murat, G.Pezzullo
// try to order routines alphabetically
// *FIXME* : need to use the assumed particle velocity instead of the speed of light
///////////////////////////////////////////////////////////////////////////////
#include "fhiclcpp/ParameterSet.h"

#include "Mu2eUtilities/inc/ModuleHistToolBase.hh"
#include "CalPatRec/inc/CalTimePeakFinder_module.hh"

// framework
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
#include "art/Utilities/make_tool.h"

// conditions
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"

#include "TrackerGeom/inc/Tracker.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"

#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"

#include "RecoDataProducts/inc/StrawHitIndex.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"

#include "Mu2eUtilities/inc/polyAtan2.hh"

using namespace std;

using CLHEP::HepVector;
using CLHEP::HepSymMatrix;
using CLHEP::Hep3Vector;

namespace mu2e {
//-----------------------------------------------------------------------------
// module constructor, parameter defaults are defiend in CalPatRec/fcl/prolog.fcl
//-----------------------------------------------------------------------------
  CalTimePeakFinder::CalTimePeakFinder(fhicl::ParameterSet const& pset) :
    art::EDFilter{pset},
    _diagLevel       (pset.get<int>            ("diagLevel"                      )),
    _debugLevel      (pset.get<int>            ("debugLevel"                     )),
    _printfreq       (pset.get<int>            ("printFrequency"                 )),
    _useAsFilter     (pset.get<int>            ("useAsFilter"                    )),
    _shLabel         (pset.get<string>         ("StrawHitCollectionLabel"        )),
    _shfLabel        (pset.get<string>         ("StrawHitFlagCollectionLabel"    )),
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
    _beta            (pset.get<double>         ("beta"                           )),  
    _dtoffset        (pset.get<double>         ("dtOffset"                       ))
  {
    consumes<ComboHitCollection>(_shLabel);
    consumes<CaloClusterCollection>(_ccmLabel);
    produces<TimeClusterCollection>();
    // produces<CalTimePeakCollection>();

    if (_debugLevel != 0) _printfreq = 1;

    if (_diagLevel  != 0) _hmanager = art::make_tool<ModuleHistToolBase>(pset.get<fhicl::ParameterSet>("diagPlugin"));
    else                  _hmanager = std::make_unique<ModuleHistToolBase>();

    _sinPitch              = sin(_pitchAngle);

    _data.minClusterEnergy =  _minClusterEnergy;
    _data.minNHits         =  _minNHits;
    _data.ntc              = 0;
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
    mu2e::GeomHandle<mu2e::Tracker> th;
    _tracker = th.get();

    mu2e::GeomHandle<mu2e::Calorimeter> ch;
    _calorimeter = ch.get();

    return true;
  }

//-----------------------------------------------------------------------------
// find input things
//-----------------------------------------------------------------------------
  bool CalTimePeakFinder::findData(const art::Event& evt) {

    auto chcolH = evt.getValidHandle<ComboHitCollection>(_shLabel);
    if (chcolH.product() != 0){
      _data.chcol = chcolH.product();
    }
    else {
      _data.chcol  = 0;
      printf(" >>> ERROR in CalTimePeakFinder::findData: ComboHitCollection with label=%s not found.\n",
             _shLabel.data());
    }

    if (evt.getByLabel(_ccmLabel, _ccH)) {
      _data.ccCollection = _ccH.product();
    }
    else {
      _data.ccCollection = 0;
      printf(" >>> ERROR in CalTimePeakFinder::findData: CaloClusterCollection with label=%s not found.\n",
             _ccmLabel.data());
    }

    return (_data.chcol != 0) && (_data.ccCollection != 0);
  }

//-----------------------------------------------------------------------------
// event entry point
//-----------------------------------------------------------------------------
  bool CalTimePeakFinder::filter(art::Event& event) {
    const char*               oname = "CalTimePeakFinder::filter";

                                        // event printout
    _iev     = event.id().event();
    if ((_debugLevel > 0) && (_iev%_printfreq) == 0) printf("[%s] : START event number %8i\n", oname,_iev);

    _data._event = &event;

    unique_ptr<TimeClusterCollection>  tcColl(new TimeClusterCollection);

    _data._tcColl = tcColl.get();

    bool ok = findData(event);

    if (ok) findTimePeaks(*_data._tcColl);
    else    printf("%s ERROR: No straw hits found in event %i\n",oname,_iev);

    int ntc = tcColl->size();

    if (_diagLevel > 0) {
//-----------------------------------------------------------------------------
// diagnostics followed by memory cleanup
//-----------------------------------------------------------------------------
      _hmanager->fillHistograms(&_data);
      _data.dtvec.clear();
    }
//-----------------------------------------------------------------------------
// put reconstructed time peaks into the event record
//-----------------------------------------------------------------------------
    event.put(std::move(tcColl));
//-----------------------------------------------------------------------------
// filtering, if requested
//-----------------------------------------------------------------------------
    if (_useAsFilter == 0) return true;

    return (ntc > 0) ;
  }
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  void CalTimePeakFinder::endJob(){ }

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  void CalTimePeakFinder::findTimePeaks(TimeClusterCollection& TimeClusterColl) {

    //    const char* oname = "CalTimePeakFinder::findTimePeaks";

    int                 ncl, nch;
    double              time, dt, tof, zstraw, cl_time;//, stime;
    double              xcl, ycl, zcl/*, dz_cl*/;
    const CaloCluster*  cl;
    const ComboHit*     hit;
    // const Straw*        straw;
    Hep3Vector          gpos, tpos;
    vector<float>       dtvec;

    //    using namespace boost::accumulators;
    static const double pi(M_PI);
    static const double twopi(2*pi);

    double              mphi(-9999.);
    double              meanDriftTime = 1.25/0.06;// half straw tube radius / drift velocity

//-----------------------------------------------------------------------------
// Loop over calorimeter clusters
//-----------------------------------------------------------------------------
    nch   = _data.chcol->size();
    ncl   = _data.ccCollection->size();

    for (int ic=0; ic<ncl; ic++) {
      cl      = &_data.ccCollection->at(ic);
      int   nsh(0);

      if ( cl->energyDep() >= _minClusterEnergy) {

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
          mphi    = polyAtan2(ycl, xcl);

          // create time peak
          TimeCluster tpeak;
//-----------------------------------------------------------------------------
// record hits in time with each peak, and accept them if they have a minimum # of hits
//-----------------------------------------------------------------------------
          for(int istr=0; istr<nch;++istr) {

            hit    = &_data.chcol->at(istr);
            time   = hit->time();
            zstraw = hit->pos().z();
//-----------------------------------------------------------------------------
// estimate time-of-flight and calculate residual between the predicted and the hit times
// 2017-03-31 P.M.: this assumes electron (e^- or e^+), not muon
//-----------------------------------------------------------------------------
            tof = (zcl-zstraw)/_sinPitch/(CLHEP::c_light*_beta);
            dt  = cl_time-(time+tof-meanDriftTime);
//--------------------------------------------------------------------------------
// check the angular distance from the calorimeter cluster
//--------------------------------------------------------------------------------
            if ((dt < _maxdt) && (dt >= _mindt)){
              double dphi = polyAtan2(hit->pos().y(), hit->pos().x()) - mphi;//phi() - mphi;

              if (dphi >  pi) dphi -= twopi;
              if (dphi < -pi) dphi += twopi;

              if (fabs(dphi) <= pi/2.) {
                tpeak._strawHitIdxs.push_back( StrawHitIndex(istr) );
                nsh += hit->nStrawHits();

		if (_diagLevel > 0) {
//-----------------------------------------------------------------------------
// accumulate diag data
//-----------------------------------------------------------------------------
		  dtvec.push_back(dt);
		}
              }
            }
          }

          if (nsh >= _data.minNHits) {
            tpeak._nsh              = nsh;
            tpeak._t0               = TrkT0(cl_time, 0.1); //dummy value for errT0
            tpeak._pos              = tpos;
            tpeak._caloCluster      = art::Ptr<mu2e::CaloCluster>(_ccH, ic);
            TimeClusterColl.push_back(tpeak);
	    _data.ntc              += 1;

	    if (_diagLevel > 0) {
	      int nh = dtvec.size();
	      for (int i=0; i<nh; i++) {
		_data.dtvec.push_back(dtvec[i]);
	      }
	    }
          }

	  if (_diagLevel > 0) dtvec.clear();
        }
      }
    }
  }

}

using mu2e::CalTimePeakFinder;
DEFINE_ART_MODULE(CalTimePeakFinder);
