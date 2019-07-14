//////////////////////////////////////////////////////////////////////////////
// framework
//
// parameter defaults: CalPatRec/fcl/prolog.fcl
//////////////////////////////////////////////////////////////////////////////
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "GeometryService/inc/GeomHandle.hh"
#include "art/Framework/Core/EDProducer.h"
#include "GeometryService/inc/DetectorSystem.hh"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
// root
#include "TVector2.h"
// data
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StereoHit.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"

// diagnostics

#include "CalPatRec/inc/DeltaFinder_types.hh"

// #include "CalPatRec/inc/LsqSums2.hh"
#include "Mu2eUtilities/inc/ModuleHistToolBase.hh"
#include "art/Utilities/make_tool.h"

#include <algorithm>
#include <cmath>
#include "CLHEP/Vector/ThreeVector.h"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"
#include "Mu2eUtilities/inc/polyAtan2.hh"

using namespace std;
using CLHEP::Hep3Vector;

namespace mu2e {

  using namespace DeltaFinderTypes;

  using CLHEP::Hep3Vector;

  class DeltaFinder: public art::EDProducer {
  public:

  protected:
    struct ChannelID {
      int Station;
      int Plane;
      int Face;
      int Panel;
      int Layer;
    };
    //-----------------------------------------------------------------------------
    // talk-to parameters: input collections and algorithm parameters
    //-----------------------------------------------------------------------------
    art::ProductToken<StrawHitCollection> const _shToken;
    art::ProductToken<ComboHitCollection> const _chToken;
    art::ProductToken<TimeClusterCollection> const _tpeakToken;
    int                                 _useTimePeaks;
    double                              _minCaloDt;
    double                              _maxCaloDt;
    double                              _pitchAngle;
    art::InputTag                       _mcDigisTag;
    float                               _minHitTime;           // min hit time
    int                                 _minNFacesWithHits;    // per station per seed
    int                                 _minNSeeds;            // min number of seeds in the delta electron cluster
    float                               _maxElectronHitEnergy; //
    float                               _minT;
    float                               _maxT;
    float                               _maxChi2Stereo;        //
    float                               _maxChi2Neighbor;      //
    float                               _maxChi2Radial;        //
    float                               _maxChi2Tot;           //
    float                               _maxDxy;
    int                                 _maxGap;
    float                               _sigmaR;
    float                               _maxDriftTime;
    float                               _maxStrawDt;
    float                               _maxDtDs;              // low-P electron travel time between two stations
    int                                 _writeStrawHits;
    int                                 _filter;

    int                                 _debugLevel;
    int                                 _diagLevel;
    int                                 _testOrder;
    std::unique_ptr<ModuleHistToolBase> _hmanager;
//-----------------------------------------------------------------------------
// cache event/geometry objects
//-----------------------------------------------------------------------------
    const ComboHitCollection*           _chcol ;
    const TimeClusterCollection*        _tpeakcol;
    StrawHitFlagCollection*             _bkgfcol;  // output collection

    const Tracker*                      _tracker;
    const DiskCalorimeter*              _calorimeter;

    float                               _tdbuff; // following Dave - time division buffer

    DeltaFinderTypes::Data_t            _data;              // all data used
    int                                 _testOrderPrinted;

    float                               _stationToCaloTOF[2][20];
    float                               _faceTOF[80];
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  public:
    explicit DeltaFinder(fhicl::ParameterSet const&);

  private:

    void         orderID  (ChannelID* X, ChannelID* Ordered);
    void         deOrderID(ChannelID* X, ChannelID* Ordered);
    void         testOrderID  ();
    void         testdeOrderID();

    bool         findData     (const art::Event&  Evt);

    int          orderHits ();

    void         findSeeds (int Station, int Face);
    void         findSeeds ();

    void         getNeighborHits(DeltaSeed* Seed, int Face1, int Face2, PanelZ_t* panelz);

    void         pruneSeeds     (int Station);

    int          checkDuplicates(int Station,
                                 int Face1, const HitData_t* Hit1,
                                 int Face2, const HitData_t* Hit2);

    void         connectSeeds      ();

    int          recoverStation    (DeltaCandidate* Delta, int Station);
    int          recoverMissingHits();

    void         runDeltaFinder();

    int          findIntersection(const HitData_t* Hit1, const HitData_t* Hit2, Intersection_t* Result);
//-----------------------------------------------------------------------------
// overloaded methods of the module class
//-----------------------------------------------------------------------------
    void beginJob() override;
    void beginRun(art::Run& ARun) override;
    void produce(art::Event& e) override;
  };

  //-----------------------------------------------------------------------------
  DeltaFinder::DeltaFinder(fhicl::ParameterSet const& pset):
    art::EDProducer{pset},
    _shToken{consumes<StrawHitCollection>(pset.get<string>("strawHitCollectionTag"))},
    _chToken{consumes<ComboHitCollection>(pset.get<string>("comboHitCollectionTag"))},
    _tpeakToken{consumes<TimeClusterCollection>(pset.get<string>("timePeakCollectionTag"))},
    _useTimePeaks          (pset.get<int>          ("useTimePeaks"                 )),
    _minCaloDt             (pset.get<double>       ("minCaloDt"                    )),
    _maxCaloDt             (pset.get<double>       ("maxCaloDt"                    )),
    _pitchAngle            (pset.get<double>       ("particleMeanPitchAngle"       )),
    _mcDigisTag            (pset.get<art::InputTag>("strawDigiMCCollectionTag"     )),

    _minHitTime            (pset.get<float>        ("minHitTime"                   )),
    _minNFacesWithHits     (pset.get<int>          ("minNFacesWithHits"            )),
    _minNSeeds             (pset.get<int>          ("minNSeeds"                    )),
    _maxElectronHitEnergy  (pset.get<float>        ("maxElectronHitEnergy"         )),
    _minT                  (pset.get<double>       ("minimumTime"                  )), // nsec
    _maxT                  (pset.get<double>       ("maximumTime"                  )), // nsec
    _maxChi2Stereo         (pset.get<float>        ("maxChi2Stereo"                )),
    _maxChi2Neighbor       (pset.get<float>        ("maxChi2Neighbor"              )),
    _maxChi2Radial         (pset.get<float>        ("maxChi2Radial"                )),
    _maxChi2Tot            (pset.get<float>        ("maxChi2Tot"                   )),
    _maxDxy                (pset.get<float>        ("maxDxy"                       )),
    _maxGap                (pset.get<int>          ("maxGap"                       )),
    _sigmaR                (pset.get<float>        ("sigmaR"                       )),
    _maxDriftTime          (pset.get<float>        ("maxDriftTime"                 )),
    _maxStrawDt            (pset.get<float>        ("maxStrawDt"                   )),
    _maxDtDs               (pset.get<float>        ("maxDtDs"                      )),
    _writeStrawHits        (pset.get<int>          ("writeStrawHits"               )),
    _filter                (pset.get<int>          ("filter"                       )),

    _debugLevel            (pset.get<int>          ("debugLevel"                   )),
    _diagLevel             (pset.get<int>          ("diagLevel"                    )),
    _testOrder             (pset.get<int>          ("testOrder"                    ))
  {
    consumesMany<ComboHitCollection>(); // Necessary because fillStrawHitIndices calls getManyByType.

    produces<StrawHitFlagCollection>("ComboHits");
    if(_writeStrawHits == 1) produces<StrawHitFlagCollection>("StrawHits");
    if (_filter) produces<ComboHitCollection>();

    _testOrderPrinted = 0;
    _tdbuff           = 80.; // mm ... abit less than 1 ns

    if (_diagLevel != 0) _hmanager = art::make_tool<ModuleHistToolBase>(pset.get<fhicl::ParameterSet>("diagPlugin"));
    else                 _hmanager = std::make_unique<ModuleHistToolBase>();
  }

  //-----------------------------------------------------------------------------
  void DeltaFinder::beginJob() {
    if (_diagLevel > 0) {
      art::ServiceHandle<art::TFileService> tfs;
      _hmanager->bookHistograms(tfs);
    }
  }

//-----------------------------------------------------------------------------
// create a Z-ordered map of the tracker
//-----------------------------------------------------------------------------
  void DeltaFinder::beginRun(art::Run& aRun) {
    mu2e::GeomHandle<mu2e::Tracker> tHandle;
    _tracker      = tHandle.get();
    _data.tracker = _tracker;

    mu2e::GeomHandle<mu2e::DiskCalorimeter> ch;
    _calorimeter = ch.get();

    ChannelID cx, co;
    int       nDisks    = _calorimeter->nDisk();
    double    disk_z[2] = {0};//given in the tracker frame

    for (int i=0; i<nDisks; ++i){
      Hep3Vector gpos = _calorimeter->disk(i).geomInfo().origin();
      Hep3Vector tpos = _calorimeter->geomUtil().mu2eToTracker(gpos);
      disk_z[i] = tpos.z();
    }

    float     z_tracker_center(0.);
    int       nPlanesPerStation(2);
    double    station_z(0);
    
    for (int planeId=0; planeId<_tracker->nPlanes(); planeId++) {
      const Plane* pln = &_tracker->getPlane(planeId);
      int  ist = planeId/nPlanesPerStation;
      int  ipl = planeId % nPlanesPerStation;
      //calculate the time-of-flight between the station and each calorimeter disk
      //for a typical Conversion Electron
      if (ipl == 0) {
	station_z = pln->origin().z();
      }else {
	station_z = (station_z + pln->origin().z())/2.;
	for (int iDisk=0; iDisk<nDisks; ++iDisk){
	  _stationToCaloTOF[iDisk][ist] = (disk_z[iDisk] - station_z)/sin(_pitchAngle)/CLHEP::c_light;
	}
      }
      
      for (int ipn=0; ipn<pln->nPanels(); ipn++) {
	const Panel* panel = &pln->getPanel(ipn);
	int face;
	if (panel->id().getPanel() % 2 == 0) face = 0;
	else                                 face = 1;
	for (int il=0; il<panel->nLayers(); ++il) {
	  cx.Station = ist;
	  cx.Plane   = ipl;
	  cx.Face    = face;
	  cx.Panel   = ipn;
	  cx.Layer   = il;
	  orderID (&cx, &co);
	  int os = co.Station; 
	  int of = co.Face;
	  int op = co.Panel;
	  PanelZ_t* pz = &_data.oTracker[os][of][op];
	  pz->fPanel = panel;
	  //-----------------------------------------------------------------------------
	  // panel caches phi of its center and the z
	  //-----------------------------------------------------------------------------
	  pz->wx  = panel->straw0Direction().x();
	  pz->wy  = panel->straw0Direction().y();
	  pz->phi = panel->straw0MidPoint().phi();
	  pz->z   = (panel->getStraw(0).getMidPoint().z()+panel->getStraw(1).getMidPoint().z())/2.;
	  int  uniqueFaceId = ipl*mu2e::StrawId::_nfaces + of;
	  _faceTOF[uniqueFaceId] = (z_tracker_center - pz->z)/sin(_pitchAngle)/CLHEP::c_light;
	}
      }	
      _data.stationUsed[ist] = 1;
      // if ((ist == 6) || (ist == 13))  _data.stationUsed[ist] = 0; // this is the past
    }
//-----------------------------------------------------------------------------
// it is enough to print that once
//-----------------------------------------------------------------------------
    if (_testOrder && (_testOrderPrinted == 0)) {
      testOrderID  ();
      testdeOrderID();
      _testOrderPrinted = 1;
    }

    if (_diagLevel != 0) _hmanager->debug(&_data,1);
  }

//------------------------------------------------------------------------------
// I'd love to use the hit flags, however that is confusing:
// - hits with very large deltaT get placed to teh middle of the wire and not flagged,
// - however, some hits within the fiducial get flagged with the ::radsel flag...
// use only "good" hits 
//-----------------------------------------------------------------------------
  int DeltaFinder::orderHits() {
    ChannelID cx, co;

    int nhits = _chcol->size();
    for (int h=0; h<nhits; ++h) {
      const ComboHit*         sh  = &(*_chcol)[h];

      if (sh->energyDep() > _maxElectronHitEnergy)         continue;
      if ( (sh->time() < _minT) || (sh->time() > _maxT) )  continue;

      cx.Station                 = sh->strawId().station();  //straw->id().getStation();
      cx.Plane                   = sh->strawId().plane() % 2;//straw->id().getPlane() % 2;
      cx.Face                    = -1;
      cx.Panel                   = sh->strawId().panel();    //straw->id().getPanel();
      cx.Layer                   = sh->strawId().layer();    //straw->id().getLayer();

					      // get Z-ordered location
      orderID(&cx, &co);

      int os       = co.Station;
      int of       = co.Face;
      int op       = co.Panel;
      int ol       = co.Layer;

      if (_useTimePeaks == 1) {
        bool               intime(false);
        int                nTPeaks  = _tpeakcol->size();
        double             hitTime  = sh->time();
        const CaloCluster* cl(nullptr);
        int                iDisk(-1);

        for (int i=0; i<nTPeaks; ++i){
          cl    = _tpeakcol->at(i).caloCluster().get();
          if (cl == nullptr) {
            printf(">>> DeltaFinder::orderHits() no CaloCluster found within the time peak %i\n", i);
            continue;
          }
          iDisk = cl->diskId();
          double    dt = cl->time() - (hitTime + _stationToCaloTOF[iDisk][os]);
          if ( (dt < _maxCaloDt) && (dt > _minCaloDt) ) {
            intime = true;
            break;
          }
        }
        if (!intime)                                    continue;
      }

      PanelZ_t* pz = &_data.oTracker[os][of][op];

      if ((os < 0) || (os >= kNStations     )) printf(" >>> ERROR: wrong station number: %i\n",os);
      if ((of < 0) || (of >= kNFaces        )) printf(" >>> ERROR: wrong face    number: %i\n",of);
      if ((op < 0) || (op >= kNPanelsPerFace)) printf(" >>> ERROR: wrong panel   number: %i\n",op);
      if ((ol < 0) || (ol >= 2              )) printf(" >>> ERROR: wrong layer   number: %i\n",ol);

      float sigw = sh->posRes(ComboHit::wire);// shp->posRes(StrawHitPosition::wire);
      pz->fHitData.push_back(HitData_t(sh,/*shp,straw,*/sigw));
    }

    return 0;
  }
  
//------------------------------------------------------------------------------
// try to recover hits of a 'Delta' candidate in a given 'Station'
// the delta candidate doesn't have hits in this station, check all hits here
// when predicting time, use the same value of Z for both layers of a given face
//-----------------------------------------------------------------------------
  int DeltaFinder::recoverStation(DeltaCandidate* Delta, int Station) {

    DeltaSeed*  new_seed (nullptr);

    for (int face=0; face<kNFaces; face++) {
      for (int ip=0; ip<kNPanelsPerFace; ip++) {
	PanelZ_t* panelz = &_data.oTracker[Station][face][ip];
	double dphi      = Delta->phi-panelz->phi;
	if (dphi < -M_PI) dphi += 2*M_PI;
	if (dphi >  M_PI) dphi -= 2*M_PI;
	if (fabs(dphi) < M_PI/3) {
//-----------------------------------------------------------------------------
// panel and seed overlap in phi, loop over hits
//-----------------------------------------------------------------------------
	  // for (int l=0; l<2; ++l) {
	    int nhits = panelz->fHitData.size();
	    for (int h=0; h<nhits; ++h) {
	      const HitData_t* hd = &panelz->fHitData[h];
	      const ComboHit*  sh = hd->fHit;
//-----------------------------------------------------------------------------
// predicted time is the partticle time, the drift time should be larger
//-----------------------------------------------------------------------------
	      if (sh->time() < Delta->T0Min(Station))                    continue;
	      if (sh->time() > Delta->T0Max(Station)+_maxDriftTime)      continue;

	      //	      if (fabs(dt) > _maxDriftTime/2 + 10)                       continue;
		  
	      double dx = sh->pos().x()-Delta->CofM.x();
	      double dy = sh->pos().y()-Delta->CofM.y();
		  
	      double dw = dx*panelz->wx+dy*panelz->wy; // distance along the wire
		  
	      double dxx = dx-panelz->wx*dw;
	      double dyy = dy-panelz->wy*dw;
		  
	      double chi2_par  = (dw*dw)/(hd->fSigW*hd->fSigW);
	      double chi2_perp = (dxx*dxx+dyy*dyy)/(_sigmaR*_sigmaR);
	      double chi2      = chi2_par + chi2_perp;
		  
	      if (chi2 >= _maxChi2Radial)          continue;
//-----------------------------------------------------------------------------
// new hit needs to be added, create a new "fake" seed for that
//-----------------------------------------------------------------------------
	      if (new_seed == NULL) new_seed = new DeltaSeed();

	      new_seed->panelz[face]  = panelz;
	      new_seed->fNHitsTot    += 1;
	      new_seed->fMaxDriftTime = _maxDriftTime;

	      if (sh->time() < new_seed->fMinTime) new_seed->fMinTime = sh->time();
	      if (sh->time() > new_seed->fMaxTime) new_seed->fMaxTime = sh->time();
	      new_seed->hitlist[face].push_back(hd);
	    }
	  // }
	}
      }
      if (new_seed) new_seed->fFaceProcessed[face] = 1;
    }
//-----------------------------------------------------------------------------
// station is processed, see if anything has been found
// some parameters of seeds found in a recovery mode are not defined because 
// there was no pre-seeding, for example
//-----------------------------------------------------------------------------
    if (new_seed) {
      new_seed->fNumber     = _data.seedHolder[Station].size();
      new_seed->fDeltaIndex = Delta->fNumber;
      _data.seedHolder[Station].push_back(new_seed);
      _data.nseeds_per_station[Station] += 1;

      Delta->seed[Station]  = new_seed;
      Delta->fNHits        += new_seed->fNHitsTot;
//-----------------------------------------------------------------------------
// update T0Min and T0Max
//-----------------------------------------------------------------------------
      if (Delta->fT0Min[Station] < new_seed->T0Min()) Delta->fT0Min[Station] = new_seed->T0Min();
      if (Delta->fT0Max[Station] > new_seed->T0Max()) Delta->fT0Max[Station] = new_seed->T0Max();

      if (Station < Delta->fFirstStation) Delta->fFirstStation = Station;
      if (Station > Delta->fLastStation ) Delta->fLastStation  = Station;
      new_seed              = nullptr;
    }

    return 0;
  }
//------------------------------------------------------------------------------
// start from looking at the "holes" in the seed pattern
// delta candidates in the list are already required to have at least 2 segments
// extend them outwards by one station
//-----------------------------------------------------------------------------
  int DeltaFinder::recoverMissingHits() {

    int ndelta = _data.deltaCandidateHolder.size();
    for (int idelta=0; idelta<ndelta; idelta++) {
      DeltaCandidate* dc = &_data.deltaCandidateHolder[idelta];
//-----------------------------------------------------------------------------
// don't extend candidates made out of one segment - but there is no such
// start from the first station to define limits
//-----------------------------------------------------------------------------
      int s1 = dc->fFirstStation;
      int s2 = dc->fLastStation-1;
      int last(-1);
      float t0min(-1.), t0max(-1.);
//-----------------------------------------------------------------------------
// first check inside "holes", skip unused stations
//-----------------------------------------------------------------------------
      for (int i=s1; i<=s2; i++) {
	if (dc->seed[i] != nullptr) {
	  last  = i; 
	  t0min = dc->fT0Min[i];
	  t0max = dc->fT0Max[i];
	  continue;
	}
	if (_data.stationUsed[i] == 0) continue;
//-----------------------------------------------------------------------------
// define expected T0 limits
//-----------------------------------------------------------------------------
	dc->fT0Min[i] = t0min-_maxDtDs*(i-last);
	dc->fT0Max[i] = t0max+_maxDtDs*(i-last);
	recoverStation(dc,i);
      }

      last  = dc->fFirstStation;
      for (int i=last-1; i>=0; i--) {
//-----------------------------------------------------------------------------
// skip empty stations
//-----------------------------------------------------------------------------
	if (_data.stationUsed[i] == 0) continue;
	dc->fT0Min[i] = dc->fT0Min[dc->fFirstStation]-_maxDtDs*(dc->fFirstStation-i);
	dc->fT0Max[i] = dc->fT0Max[dc->fFirstStation]+_maxDtDs*(dc->fFirstStation-i);
	recoverStation(dc,i);
//-----------------------------------------------------------------------------
// so far, do not allow holes while extending 
//-----------------------------------------------------------------------------
	if (dc->fFirstStation != i) break;
      }

      last = dc->fLastStation;
      for (int i=last+1; i<kNStations; i++) {
//-----------------------------------------------------------------------------
// skip empty stations
//-----------------------------------------------------------------------------
	if (_data.stationUsed[i] == 0) continue;
	dc->fT0Min[i] = dc->fT0Min[dc->fLastStation]-_maxDtDs*(i-dc->fLastStation);
	dc->fT0Max[i] = dc->fT0Max[dc->fLastStation]+_maxDtDs*(i-dc->fLastStation);
	recoverStation(dc,i);
//-----------------------------------------------------------------------------
// so far, do not allow holes while extending 
//-----------------------------------------------------------------------------
	if (dc->fLastStation != i) break;
      }
    }

    return 0;
  }

//-----------------------------------------------------------------------------
  void  DeltaFinder::runDeltaFinder() {

    for (int s=0; s<kNStations; ++s) {
      for (int f=0; f<kNFaces; ++f) {
        for (int p=0; p<3; ++p) {
          PanelZ_t* panelz = &_data.oTracker[s][f][p];
          panelz->fHitData.clear() ;
        }
      }
    }
    orderHits();
    findSeeds();
    connectSeeds();
    recoverMissingHits();
  }

//-----------------------------------------------------------------------------
  bool DeltaFinder::findData(const art::Event& Evt) {
    _chcol    = nullptr;
    _tpeakcol = nullptr;

    if (_useTimePeaks == 1){
      auto tpeakH    = Evt.getValidHandle(_tpeakToken);
      _tpeakcol      = tpeakH.product();
      _data.tpeakcol = _tpeakcol;  // FIXME
    }

    auto shH    = Evt.getValidHandle(_chToken);
    _chcol      = shH.product();
    _data.chcol = _chcol;  // FIXME

    return (_chcol != 0);
  }

//-----------------------------------------------------------------------------
  void DeltaFinder::produce(art::Event& Event) {

    if (_debugLevel) printf(">>> DeltaFinder::produce  event number: %10i\n",Event.event());  
//-----------------------------------------------------------------------------
// clear memory in the beginning of event processing and cache event pointer
//-----------------------------------------------------------------------------
    _data.event  = &Event;

    _data.nseeds = 0;
    _data.debugLevel = _debugLevel;

    for (int is=0; is<kNStations; is++) {
      _data.nseeds_per_station[is] = 0;

      for (auto ds=_data.seedHolder[is].begin(); ds!=_data.seedHolder[is].end(); ds++) {
        delete *ds;
      }
      _data.seedHolder[is].clear();
    }

    _data.deltaCandidateHolder.clear();
//-----------------------------------------------------------------------------
// process event
//-----------------------------------------------------------------------------
    if (! findData(Event)) {
      const char* message = "mu2e::DeltaFinder_module::produce: data missing or incomplete";
      throw cet::exception("RECO")<< message << endl;
    }

    runDeltaFinder();
//-----------------------------------------------------------------------------
// form output - copy input flag collection - do we need it ?
//-----------------------------------------------------------------------------
    unique_ptr<StrawHitFlagCollection> bkgfcol(new StrawHitFlagCollection());
    _bkgfcol = bkgfcol.get();
    int nch = _chcol->size();
    _bkgfcol->reserve(nch); // add initialization loop *FIXME* 
    for (int i=0; i<nch; i++) {
      StrawHitFlag flag;
      const ComboHit* sh  = &(*_chcol)[i];
      if (sh->energyDep() < _maxElectronHitEnergy) flag.merge(StrawHitFlag::energysel);
      _bkgfcol->push_back(flag);
    }
    _data.shfcol = _bkgfcol;

    const ComboHit* sh0(0);
    if (nch > 0) sh0 = &_chcol->at(0);

    StrawHitFlag deltamask(StrawHitFlag::bkg);

    int ndeltas = _data.deltaCandidateHolder.size();

    for (int i=0; i<ndeltas; i++) {
      DeltaCandidate* dc = &_data.deltaCandidateHolder.at(i);
      for (int station=dc->fFirstStation; station<=dc->fLastStation; station++) {
	DeltaSeed* ds = dc->seed[station];
	if (ds != nullptr) {
//-----------------------------------------------------------------------------
// loop over the hits and flag each of them as delta
//-----------------------------------------------------------------------------
	  for (int face=0; face<kNFaces; face++) {
	    int nh = ds->NHits(face);
	    for (int ih=0; ih<nh; ih++) {
	      const ComboHit* sh = ds->HitData(face,ih)->fHit;
	      int loc = sh-sh0;
	      _bkgfcol->at(loc).merge(deltamask);
	    }
	  }
	}
      }
    }
//-----------------------------------------------------------------------------
// in the end of event processing fill diagnostic histograms
//-----------------------------------------------------------------------------
    if (_diagLevel  > 0) _hmanager->fillHistograms(&_data);
    if (_debugLevel > 0) _hmanager->debug(&_data,2);


//-----------------------------------------------------------------------------
// create the collection of StrawHitFlag for the StrawHitCollection
//-----------------------------------------------------------------------------
    if(_writeStrawHits == 1){
      auto shH = Event.getValidHandle<StrawHitCollection>(_shToken);
      const StrawHitCollection* shcol = shH.product();
      // first, copy over the original flags
      unsigned nsh = shcol->size();
      std::unique_ptr<StrawHitFlagCollection> shfcol(new StrawHitFlagCollection(nsh));
      std::vector<std::vector<StrawHitIndex> > shids;
      _chcol->fillStrawHitIndices(Event,shids);
      for(size_t ich = 0;ich < _chcol->size();++ich) {
    	StrawHitFlag flag = bkgfcol->at(ich);
    	flag.merge((*_chcol)[ich].flag());
    	for(auto ish : shids[ich])
    	  (*shfcol)[ish] = flag;
      }

      Event.put(std::move(shfcol),"StrawHits");
    }
//-----------------------------------------------------------------------------
// create the collection of StrawHitFlag for the StrawHitCollection
//-----------------------------------------------------------------------------
    if (_filter == 1){
      auto chcol = std::make_unique<ComboHitCollection>();
      chcol->reserve(nch);
      // same parent as the original collection
      chcol->setParent(_chcol->parent());
      for(int ich=0;ich < nch; ++ich){
        StrawHitFlag const& flag = bkgfcol->at(ich);
        if(!flag.hasAnyProperty(StrawHitFlag::bkg)) {
          chcol->push_back((*_chcol)[ich]);
          chcol->back()._flag.merge(flag);
        }
      }
      Event.put(std::move(chcol));
    }
//-----------------------------------------------------------------------------
// finally, put the output flag collection into the event
//-----------------------------------------------------------------------------
    Event.put(std::move(bkgfcol),"ComboHits");
  }

//-----------------------------------------------------------------------------
  void DeltaFinder::orderID(ChannelID* X, ChannelID* O) {
    if (X->Panel % 2 == 0) X->Face = 0;
    else                   X->Face = 1; // define original face

    O->Station = X->Station; // stations already ordered
    O->Plane   = X->Plane;   // planes already ordered, but not necessary for ordered construct

    if (X->Station % 2 == 0) {
      if (X->Plane == 0) O->Face = 1 - X->Face;
      else               O->Face = X->Face + 2;
    }
    else {
      if (X->Plane == 0) O->Face = X->Face;
      else               O->Face = 3 - X->Face; // order face
    }

    O->Panel = int(X->Panel/2);                // order panel

    int n = X->Station + X->Plane + X->Face;   // pattern has no intrinsic meaning, just works
    if (n % 2 == 0) O->Layer = 1 - X->Layer;
    else            O->Layer = X->Layer;       // order layer
  }
  
//-----------------------------------------------------------------------------
  void DeltaFinder::deOrderID(ChannelID* X, ChannelID* O) {

    X->Station = O->Station;

    X->Plane   = O->Plane;

    if(O->Station % 2 ==  0) {
      if(O->Plane == 0) X->Face = 1 - O->Face;
      else X->Face = O->Face - 2;
    }
    else {
      if(O->Plane == 0) X->Face = O->Face;
      else X->Face = 3 - O->Face;
    }

    if(X->Face == 0) X->Panel = O->Panel * 2;
    else X->Panel = 1 + (O->Panel * 2);

    int n = X->Station + X->Plane + X->Face;
    if(n % 2 == 0) X->Layer = 1 - O->Layer;
    else X->Layer = O->Layer;
  }
 
//-----------------------------------------------------------------------------
// testOrderID & testdeOrderID not used in module, only were used to make sure OrderID and deOrderID worked as intended   
//-----------------------------------------------------------------------------
  void DeltaFinder::testOrderID() {

    ChannelID x, o;

    for (int s=0; s<2; ++s) {
      for (int pl=0; pl<2; ++pl) {
	for (int pa=0; pa<6; ++pa) {
	  for (int l=0; l<2; ++l) {
	    x.Station = s;
	    x.Plane   = pl;
	    x.Panel   = pa;
	    x.Layer   = l;
	    orderID(&x, &o);
	    printf(" testOrderID: Initial(station = %i, plane = %i, face = %i, panel = %i, layer = %i)", 
		   x.Station, x.Plane, x.Face, x.Panel, x.Layer);
	    printf("  Ordered(station = %i, plane = %i, face = %i, panel = %i, layer = %i)\n",
		   o.Station, o.Plane, o.Face, o.Panel, o.Layer);
	  }
	}
      }
    }
  }

//-----------------------------------------------------------------------------
  void DeltaFinder::testdeOrderID() {

    ChannelID x, o;

    for (int s=0; s<2; ++s) {
      for (int f=0; f<4; ++f) {
        for (int pa=0; pa<3; ++pa) {
          for (int l=0; l<2; ++l) {

            o.Station          = s;
            o.Face             = f;
            if (f < 2) o.Plane = 0;
            else       o.Plane = 1;
            o.Panel            = pa;
            o.Layer            = l;

            deOrderID(&x, &o);

            printf(" testdeOrderID: Initial(station = %i, plane = %i, face = %i, panel = %i, layer = %i)",
                   x.Station,x.Plane,x.Face,x.Panel,x.Layer);
            printf("  Ordered(station = %i, plane = %i, face = %i, panel = %i, layer = %i)\n",
                   o.Station,o.Plane,o.Face,o.Panel,o.Layer);
          }
        }
      }
    }
  }

//-----------------------------------------------------------------------------
// make sure the two hits used to make a new seed are not a part of an already found seed
//-----------------------------------------------------------------------------
  int DeltaFinder::checkDuplicates(int Station, int Face1, const HitData_t* Hit1, int Face2, const HitData_t* Hit2) {

    bool h1_found(false), h2_found(false);

    int nseeds = _data.seedHolder[Station].size();
    for (int i=0; i<nseeds; i++) {
      DeltaSeed* seed = _data.seedHolder[Station][i];

      int nhits = seed->hitlist[Face1].size();
      for (int ih=0; ih<nhits; ih++) {
        const HitData_t* hit = seed->hitlist[Face1][ih];
        if (hit == Hit1) {
          h1_found = true;
          break;
        }
      }

      if (! h1_found) continue;  // with the next seed
//-----------------------------------------------------------------------------
// Hit1 was found, search the same seed for Hit2
//-----------------------------------------------------------------------------
      nhits = seed->hitlist[Face2].size();
      for (int ih=0; ih<nhits; ih++) {
        const HitData_t* hit = seed->hitlist[Face2][ih];
        if (hit == Hit2) {
          h2_found = true;
          break;
        }
      }

      if (h2_found) {
//-----------------------------------------------------------------------------
// both Hit1 and Hit2 were found within the same seed
//-----------------------------------------------------------------------------
	return 1;
      }
    }
    return 0;
  }


//-----------------------------------------------------------------------------
  int DeltaFinder::findIntersection(const HitData_t* Hd1, const HitData_t* Hd2, Intersection_t* Result) {
    double x1, y1, x2, y2, nx1, ny1, nx2, ny2;

    const Hep3Vector& p1 = Hd1->fHit->centerPosCLHEP();//fStraw->getMidPoint();

    x1 =  p1.x();
    y1 =  p1.y();

    const Hep3Vector& p2 = Hd2->fHit->centerPosCLHEP();//fStraw->getMidPoint();
    x2 =  p2.x();
    y2 =  p2.y();

    const Hep3Vector& wdir1 = Hd1->fHit->centerPosCLHEP();//fStraw->getDirection();
    nx1 = wdir1.x();
    ny1 = wdir1.y();

    const Hep3Vector& wdir2 = Hd2->fHit->centerPosCLHEP();//fStraw->getDirection();
    nx2 = wdir2.x();
    ny2 = wdir2.y();

    double n1n2  = nx1*nx2+ny1*ny2;
    double r12n1 = (x1-x2)*nx1+(y1-y2)*ny1;
    double r12n2 = (x1-x2)*nx2+(y1-y2)*ny2;
//-----------------------------------------------------------------------------
// t1 and t2 are distances to the intersection point from the centers of the 
// corresponding wires
//-----------------------------------------------------------------------------
    Result->t1 = (n1n2*r12n2-r12n1)/(1-n1n2*n1n2);
    Result->t2 = (r12n2-n1n2*r12n1)/(1-n1n2*n1n2);

					// in 2D, the lines intersect, take one
    Result->x = x1+nx1*Result->t1;
    Result->y = y1+ny1*Result->t1;
//-----------------------------------------------------------------------------
// now define distances to the hits
//-----------------------------------------------------------------------------
    Hep3Vector h1 = Hd1->fHit->posCLHEP();
    Result->wd1 = (h1.x()-Result->x)*nx1+(h1.y()-Result->y)*ny1;
    Hep3Vector h2 = Hd2->fHit->posCLHEP();
    Result->wd2 = (h2.x()-Result->x)*nx2+(h2.y()-Result->y)*ny2;

    return 0;
  }


//-----------------------------------------------------------------------------
// pick up neighboring hits in 'Face'
//-----------------------------------------------------------------------------
  void DeltaFinder::getNeighborHits(DeltaSeed* Seed, int Face, int Face2, PanelZ_t* panelz) {

    DeltaFinderTypes::Intersection_t     res;

    vector<HitData_t*> hits;

    const HitData_t* hd1 = Seed->HitData(Face,0);
    // const Straw* straw1  = hd1->fStraw;
    float minrad         = hd1->fHit->centerPosCLHEP().perp();//straw1->getMidPoint().perp();
    float maxrad         = minrad;

    // for (int l=0; l<2; ++l) {
    int nh = panelz->fHitData.size();
    for (int h=0; h<nh; ++h) {
      HitData_t* hd  = &panelz->fHitData[h];
      if (hd == hd1) continue ;
      const ComboHit* sh = hd->fHit;

      if (sh->time()-Seed->T0Max() > _maxDriftTime          ) continue;
      if (sh->time()               < Seed->T0Min()          ) continue;

      hd->fDr        = fabs(hd->fRMid-minrad);
      hits.push_back(hd);
    }
    // }
    //-----------------------------------------------------------------------------
    // sorted hits in dr wrt the seed hit.
    // I know that sorting could be done with more elegance, don't care at the moment
    // also the diagnostics data could be fully separated from the data itself
    //-----------------------------------------------------------------------------
    int nhits = hits.size();
    for (int i=0; i<nhits-1; i++) {
      HitData_t** hi = &hits[i];
      for (int j=i+1; j<nhits; j++) {
        HitData_t** hj = &hits[j];
        if ((*hi)->fDr >= (*hj)->fDr) {
          HitData_t* h = hits[i];
          hits[i]      = hits[j];
          hits[j]      = h;
        }
      }
    }
//-----------------------------------------------------------------------------
// hits are orders in distance from the first pre-seed hit
// loop over them again, 'hd2' - the second seed hit
//-----------------------------------------------------------------------------
    const HitData_t* hd2 = Seed->HitData(Face2,0);

    for (int i=0; i<nhits; i++) {
      const HitData_t* hd = hits[i];
      float rad     = hd->fRMid;
      if ((minrad-rad > 10) || (rad-maxrad > 10) ) continue;
//-----------------------------------------------------------------------------
// radially we're OK, check distance from the intersection
//-----------------------------------------------------------------------------
      CLHEP::Hep3Vector pos;

      DeltaFinderTypes::findIntersection(hd, hd2, &res);

      float chi = res.wd1/hd->fSigW;

      if (chi*chi < _maxChi2Neighbor) {
//-----------------------------------------------------------------------------
// OK along the wire, add hit as a neighbor
//-----------------------------------------------------------------------------
	if (rad < minrad) minrad = rad;
	else              maxrad = rad;

	Seed->hitlist[Face].push_back(hd);

	const ComboHit* sh = hd->fHit;
	if (sh->time() < Seed->fMinTime) Seed->fMinTime = sh->time();
	if (sh->time() > Seed->fMaxTime) Seed->fMaxTime = sh->time();
      }
    }
  }


//-----------------------------------------------------------------------------
// find delta electron seeds in 'Station' with hits in faces 'f' and 'f+1'
// do not consider proton hits with eDep > _minHitEnergy
//-----------------------------------------------------------------------------
  void DeltaFinder::findSeeds(int Station, int Face) {

    DeltaFinderTypes::Intersection_t  res;

    for (int p=0; p<3; ++p) {                        // loop over panels
      PanelZ_t* panelz = &_data.oTracker[Station][Face][p];
      // for (int l=0; l<2; ++l) {                     // loop over layers
      int hitsize1 = panelz->fHitData.size();
      for (int h1=0; h1<hitsize1; ++h1) {             // loop over hits/hit positions
        //-----------------------------------------------------------------------------
        // hit has not been used yet to start a seed,
        // however it could've been used as a second seed
        //-----------------------------------------------------------------------------
        HitData_t* hd1 = &panelz->fHitData[h1];
        const ComboHit* sh = hd1->fHit;
        //        if (sh->energyDep() >= _maxElectronHitEnergy)  continue;
        // if (fabs(sh->dt())  >= _maxStrawDt          )  continue;//FIXME!
        float ct = sh->time();
        if (ct              <  _minHitTime          )  continue;
        // const Straw* straw1 = hd1->fStraw;
        int counter         = 0;                // number of stereo candidates hits close to set up counter
        //-----------------------------------------------------------------------------
        // loop over the second faces
        //-----------------------------------------------------------------------------
        for (int f2=Face+1; f2<kNFaces; f2++) {
          for (int p2=0; p2<3; ++p2) {         // loop over panels
            PanelZ_t* panelz2 = &_data.oTracker[Station][f2][p2];
            //-----------------------------------------------------------------------------
            // check if the two panels overlap in XY
            // 2D angle between the vectors pointing to the panel centers, can't be greater than pi
            //-----------------------------------------------------------------------------
            float dphi = panelz2->phi - panelz->phi;
            if (dphi < -M_PI) dphi += 2*M_PI;
            if (dphi >  M_PI) dphi -= 2*M_PI;
            if (abs(dphi) >= 2*M_PI/3.) continue;
            //-----------------------------------------------------------------------------
            // panels do overlap
            //-----------------------------------------------------------------------------
            // for (int l2=0; l2<2;++l2) {
            int hitsize2 = panelz2->fHitData.size();
            for (int h2=0; h2<hitsize2;++h2) {
              HitData_t* hd2 = &panelz2->fHitData[h2];
              const ComboHit* sh2 = hd2->fHit;
              //                  if (sh2->energyDep() >= _maxElectronHitEnergy)  continue;
              // if (fabs(sh2->dt())  >= _maxStrawDt)            continue;//FIXME!
              float t2 = sh2->time();
              if (t2 < _minHitTime)                           continue;
              float dt = abs(t2 - ct);
              if (dt >= _maxDriftTime)                        continue;
              ++counter;                                            // number of hits close to the first one
              //-----------------------------------------------------------------------------
              // intersect the two straws, we need coordinates of the intersection point and
              // two distances from hits to the intersection point, 4 numbers in total
              //-----------------------------------------------------------------------------
              DeltaFinderTypes::findIntersection(hd1,hd2,&res);
              //-----------------------------------------------------------------------------
              // make sure the two straws do intersect
              //-----------------------------------------------------------------------------
              // if (fabs(res.t1) >= straw1->getHalfLength())    continue;

              // const Straw* straw2              = hd2->fStraw;
              // if (fabs(res.t2) >= straw2->getHalfLength())    continue;
              //-----------------------------------------------------------------------------
              // both StrawHit's are required to be close enough to the intersection point
              //-----------------------------------------------------------------------------
              float chi1 = res.wd1/hd1->fSigW;
              if (chi1*chi1 >= _maxChi2Stereo)                continue;
              float chi2 = res.wd2/hd2->fSigW;
              if (chi2*chi2 >= _maxChi2Stereo)                continue;
              //-----------------------------------------------------------------------------
              // check whether there already is a seed containing both hits
              //-----------------------------------------------------------------------------
              int is_duplicate = checkDuplicates(Station,Face,hd1,f2,hd2);
              if (is_duplicate)                               continue;
              //-----------------------------------------------------------------------------
              // new seed
              //-----------------------------------------------------------------------------
              DeltaSeed* seed = new DeltaSeed();
              seed->fStation             =  Station;
              seed->fNumber              =  _data.seedHolder[Station].size();
              seed->fType                = 10*Face+f2;
              seed->fNFacesWithHits      = 2;
              seed->fFaceProcessed[Face] = 1;
              seed->fFaceProcessed[f2  ] = 1;
              seed->fMaxDriftTime        = _maxDriftTime;
              // could these be redefined? - in principle, yes...
              hd1->fChi2Min    = chi1*chi1;
              hd2->fChi2Min    = chi2*chi2;
              hd1->fSeedNumber = seed->fNumber;
              hd2->fSeedNumber = seed->fNumber;

              seed->fMinTime = sh->time();
              if (sh2->time() > seed->fMinTime) {
                seed->fMaxTime = sh2->time();
              }
              else {
                seed->fMinTime = sh2->time();
                seed->fMaxTime = sh->time();
              }

              seed->hitlist[Face].push_back(hd1);
              seed->hitlist[f2  ].push_back(hd2);

              // getNeighborHits(seed, Face, f2, panelz);

              // CLHEP::Hep3Vector smpholder(straw1->getMidPoint()); // use precalculated
              CLHEP::Hep3Vector         smpholder(hd1->fHit->centerPosCLHEP());

              int nh1 = seed->hitlist[Face].size();
              for(int h3=1; h3<nh1; ++h3) {
                smpholder  += seed->hitlist[Face][h3]->fHit->centerPosCLHEP();//fStraw->getMidPoint();
              }

              // getNeighborHits(seed, f2, Face, panelz2);

              // CLHEP::Hep3Vector smp2holder(straw2->getMidPoint());
              CLHEP::Hep3Vector          smp2holder(hd2->fHit->centerPosCLHEP());
              int nh2 = seed->hitlist[f2].size();
              for(int h4=1; h4<nh2; ++h4) {
                smp2holder += seed->hitlist[f2][h4]->fHit->centerPosCLHEP();//fStraw->getMidPoint();
              }

              CLHEP::Hep3Vector CofMsmp1 = smpholder /nh1;
              CLHEP::Hep3Vector CofMsmp2 = smp2holder/nh2;

              const Hep3Vector& CofMdir1 = hd1->fHit->wdirCLHEP();//straw1->getDirection();
              const Hep3Vector& CofMdir2 = hd2->fHit->wdirCLHEP();//straw2->getDirection();
              TwoLinePCA pca(CofMsmp1, CofMdir1, CofMsmp2, CofMdir2);

              seed->CofM         = 0.5*(pca.point1() + pca.point2());
              seed->fHitData[0]  = hd1;
              seed->fHitData[1]  = hd2;
              seed->fChi21       = chi1*chi1;
              seed->fChi22       = chi2*chi2;
              seed->panelz[Face] = panelz;
              seed->panelz[f2]   = panelz2;
              seed->fNHitsTot    = nh1+nh2;

              _data.seedHolder[Station].push_back(seed);
              //-----------------------------------------------------------------------------
              // book-keeping: increment total number of found seeds
              //-----------------------------------------------------------------------------
              _data.nseeds                      += 1;
              _data.nseeds_per_station[Station] += 1;
            }
            // }
          }
        }
        //-----------------------------------------------------------------------------
        // this is needed for diagnostics only
        //-----------------------------------------------------------------------------
        if (_diagLevel > 0) {
          hd1->fNSecondHits  = counter ;
        }
      }
    }
    // }
  }

//-----------------------------------------------------------------------------
// some of found seeds could be duplicates or ghosts
// in case two DeltaSeeds share the first seed hit, leave only the best one
// all seeds we're loooping over have been reconstructed within the same
// also reject seeds with Chi2Tot > 10
//-----------------------------------------------------------------------------
  void DeltaFinder::pruneSeeds(int Station) {
    int nseeds =  _data.seedHolder[Station].size();

    for (int i1=0; i1<nseeds-1; i1++) {
      DeltaSeed* ds1 = _data.seedHolder[Station][i1];
      if (ds1->fGood < 0) continue;

      if (ds1->Chi2Tot() > _maxChi2Tot) {
        ds1->fGood = -1000-i1;
        continue;
      }

      float tmean1 = (ds1->fMinTime+ds1->fMaxTime)/2.;

      for (int i2=i1+1; i2<nseeds; i2++) {
	DeltaSeed* ds2 = _data.seedHolder[Station][i2];
	if (ds2->fGood < 0) continue;

	if (ds2->Chi2Tot() > _maxChi2Tot) {
	  ds2->fGood = -1000-i2;
	  continue;
	}

	float tmean2 = (ds2->fMinTime+ds2->fMaxTime)/2.;

	if (fabs(tmean1-tmean2) > _maxDriftTime) continue;
//-----------------------------------------------------------------------------
// the two segments are close in time and space, check hit overlap
// *FIXME* didn't check distance !!!!!
// so far, allow duplicates during the search
// the two DeltaSeeds share could have significantly overlapping hit content
//-----------------------------------------------------------------------------
	int noverlap = 0;
	int nfaces_with_overlap = 0;
	for (int face=0; face<kNFaces; face++) {
	  int nov = 0;
	  int nh1 = ds1->hitlist[face].size();
	  for (int ih1=0; ih1<nh1; ih1++) {
	    const HitData_t* hh1 = ds1->hitlist[face][ih1];
	    int nh2 = ds2->hitlist[face].size();
	    for (int ih2=0; ih2<nh2; ih2++) {
	      const HitData_t* hh2 = ds2->hitlist[face][ih2];
	      if (hh2 == hh1) {
		nov += 1;
		break;
	      }
	    }
	  }
	  noverlap += nov;
	  if (nov != 0) nfaces_with_overlap += 1;
	}
//-----------------------------------------------------------------------------
// special treatment of 2-hit seeds to reduce the number of ghosts
//-----------------------------------------------------------------------------
	if (ds1->fNHitsTot == 2) {
	  if (nfaces_with_overlap > 0) {
	    if (ds2->fNFacesWithHits > 2) {
	      ds1->fGood = -1000-i2;
	      break;
	    }
	    else {
					// the second one also has 2 faces with hits

	      if (ds1->Chi2AllDof() <  ds2->Chi2AllDof()) ds2->fGood = -1000-i1;
	      else {
		ds1->fGood = -1000-i2;
		break;
	      }
	    }
	  }
	}

	if (ds2->fNHitsTot == 2) {
	  if (nfaces_with_overlap > 0) {
//-----------------------------------------------------------------------------
// the 2nd seed has only 2 hits and there is an overlap
//-----------------------------------------------------------------------------
	    if (ds1->fNFacesWithHits > 2)                 ds2->fGood = -1000-i1;
	    else {
					// the second one also has 2 faces with hits

	      if (ds1->Chi2AllDof() <  ds2->Chi2AllDof()) ds2->fGood = -1000-i1;
	      else {
		ds1->fGood = -1000-i2;
		break;
	      }
	    }
	  }
	}

	if (nfaces_with_overlap > 1) {
	  if ((noverlap >= ds1->fNHitsTot*0.6) || (noverlap >= 0.6*ds2->fNHitsTot)) {
//-----------------------------------------------------------------------------
// overlap significant, leave in only one DeltaSeed - which one? 
//-----------------------------------------------------------------------------
	    if      (ds1->Chi2AllDof() <  ds2->Chi2AllDof()) ds2->fGood = -1000-i1;
	    else if (ds1->Chi2AllDof() >= ds2->Chi2AllDof()) {
	      ds1->fGood = -1000-i2;
	      break;
	    }
	    else {
//-----------------------------------------------------------------------------
// the same number of hits, choose candidate with lower chi2
// should not be getting here
//-----------------------------------------------------------------------------
	      if (ds1->Chi2AllDof() < ds2->Chi2AllDof()) ds2->fGood = -1000-i1;
	      else {
		ds1->fGood = -1000-i2;
		break;
	      }
	    }
	  }
	}
      }
    }
  }
     
//-----------------------------------------------------------------------------
// TODO: update the time as more hits are added
//-----------------------------------------------------------------------------
  void DeltaFinder::findSeeds() {

    for (int s=0; s<kNStations; ++s) { 
      for (int f1=0; f1<kNFaces-1; ++f1) {
//-----------------------------------------------------------------------------
// 'last' - number of seeds found so far
//-----------------------------------------------------------------------------
	int last = _data.seedHolder[s].size();
	
	findSeeds(s,f1);
//-----------------------------------------------------------------------------
// for seeds with hits in faces (f,f+1), (f,f+2), (f,f+3) find hits in other two faces
//-----------------------------------------------------------------------------
	int nseeds = _data.seedHolder[s].size();
	for (int iseed=last; iseed<nseeds; iseed++) {
	  DeltaSeed* seed = _data.seedHolder[s][iseed];
	  double seed_phi = polyAtan2(seed->CofM.y(), seed->CofM.x());//seed->CofM.phi();              // check to find right panel
//-----------------------------------------------------------------------------
// simultaneously update CoM coordinates
//-----------------------------------------------------------------------------
	  double sx(0), sy(0), snx2(0),snxny(0), sny2(0), snxnr(0), snynr(0);

	  for (int face=f1; face<kNFaces; face++) {
	    int nh = seed->NHits(face);
	    for (int ih=0; ih<nh; ih++) {
	      const HitData_t* hd = seed->HitData(face,ih);
	      // const Straw*     s  = hd->fStraw;

	      double x0 = hd->fHit->pos().x();// CHECK IT! s->getMidPoint().x();
	      double y0 = hd->fHit->pos().y();// CHECK IT! s->getMidPoint().y();
	      double nx = hd->fHit->wdir().x();//          s->getDirection().x();
	      double ny = hd->fHit->wdir().y();//          s->getDirection().y();
	      double nr = nx*x0+ny*y0;
	      
	      sx    += x0;
	      sy    += y0;
	      snx2  += nx*nx;
	      snxny += nx*ny;
	      sny2  += ny*ny;
	      snxnr += nx*nr;
	      snynr += ny*nr;
	    }
	  }
//-----------------------------------------------------------------------------
// loop over remaining two faces, 'f2' - face in question
//-----------------------------------------------------------------------------
	  for (int f2=0; f2<kNFaces; f2++) {
	    if (seed->fFaceProcessed[f2] == 1)                              continue;
//-----------------------------------------------------------------------------
// face is different from the two first faces used
//-----------------------------------------------------------------------------
	    for (int p2=0; p2<3; ++p2) {
	      PanelZ_t* panelz = &_data.oTracker[s][f2][p2];
	      double dphi      = seed_phi-panelz->phi;
	      if (dphi < -M_PI) dphi += 2*M_PI;
	      if (dphi >  M_PI) dphi -= 2*M_PI;
	      if (fabs(dphi) >= M_PI/3)                                     continue;
//-----------------------------------------------------------------------------
// panel overlaps with the seed, look at its hits
//-----------------------------------------------------------------------------
	      // for(int l=0; l<2; ++l) {
		int psize = panelz->fHitData.size();
		for (int h=0; h<psize; ++h) { // find hit
//-----------------------------------------------------------------------------
// 2017-10-05 PM: consider all hits 
// hit time should be consistent with the already existing times - the difference
// between any two measured hit times should not exceed _maxDriftTime 
// (_maxDriftTime represents the maximal drift time in the straw, should there be some tolerance?)
//-----------------------------------------------------------------------------
		  HitData_t* hd      = &panelz->fHitData[h];
		  const ComboHit* sh = hd->fHit;

		  if (sh->time()-seed->T0Max() > _maxDriftTime          ) continue;
		  if (sh->time()               < seed->T0Min()          ) continue;

		  // const StrawHitPosition* shp  = hd->fPos;
		  CLHEP::Hep3Vector       dxyz = sh->posCLHEP()-seed->CofM;// shp->posCLHEP()-seed->CofM; // distance from hit to preseed
//-----------------------------------------------------------------------------
// split into wire parallel and perpendicular components
//-----------------------------------------------------------------------------
		  const CLHEP::Hep3Vector& wdir = hd->fHit->wdirCLHEP();//fStraw->getDirection();
		  CLHEP::Hep3Vector d_par    = (dxyz.dot(wdir))/(wdir.dot(wdir))*wdir; 
		  CLHEP::Hep3Vector d_perp_z = dxyz-d_par;
		  float  d_perp              = d_perp_z.perp();
		  double sigw                = hd->fSigW;
		  float  chi2_par            = (d_par.mag()/sigw)*(d_par.mag()/sigw);
		  float  chi2_perp           = (d_perp/_sigmaR)*(d_perp/_sigmaR);
		  float  chi2                = chi2_par + chi2_perp;
		  if (chi2 >= _maxChi2Radial)                             continue;
//-----------------------------------------------------------------------------
// add hit
//-----------------------------------------------------------------------------
		  hd->fChi2Min = chi2;
		  seed->hitlist[f2].push_back(hd);

		  if (sh->time() < seed->fMinTime) seed->fMinTime = sh->time();
		  if (sh->time() > seed->fMaxTime) seed->fMaxTime = sh->time();

		  seed->fNHitsTot++;
//-----------------------------------------------------------------------------
// in parallel, update coordinate sums
//-----------------------------------------------------------------------------
		  // const Straw* straw  = hd->fStraw;

		  double x0 = hd->fHit->pos().x();//straw->getMidPoint().x();
		  double y0 = hd->fHit->pos().y();// straw->getMidPoint().y();
		  double nx = hd->fHit->wdir().x();// straw->getDirection().x();
		  double ny = hd->fHit->wdir().y();//  straw->getDirection().y();
		  double nr = nx*x0+ny*y0;
		      
		  sx    += x0;
		  sy    += y0;
		  snx2  += nx*nx;
		  snxny += nx*ny;
		  sny2  += ny*ny;
		  snxnr += nx*nr;
		  snynr += ny*nr;
		}
	      // }
	    }
//-----------------------------------------------------------------------------
// update seed time and X and Y coordinates, accurate knowledge of Z is not very relevant
//-----------------------------------------------------------------------------
	    double x_mean, y_mean, nxny_mean, nx2_mean, ny2_mean, nxnr_mean, nynr_mean;

	    x_mean    = sx   /seed->fNHitsTot;
	    y_mean    = sy   /seed->fNHitsTot;
	    nxny_mean = snxny/seed->fNHitsTot;
	    nx2_mean  = snx2 /seed->fNHitsTot;
	    ny2_mean  = sny2 /seed->fNHitsTot;
	    nxnr_mean = snxnr/seed->fNHitsTot;
	    nynr_mean = snynr/seed->fNHitsTot;

	    double d = (1-nx2_mean)*(1-ny2_mean)-nxny_mean*nxny_mean;
	    
	    double x0 = ((x_mean-nxnr_mean)*(1-ny2_mean)+(y_mean-nynr_mean)*nxny_mean)/d;
	    double y0 = ((y_mean-nynr_mean)*(1-nx2_mean)+(x_mean-nxnr_mean)*nxny_mean)/d;

	    seed->CofM.setX(x0);
	    seed->CofM.setY(y0);

	    if (seed->hitlist[f2].size() > 0) seed->fNFacesWithHits++;
	    seed->fFaceProcessed[f2] = 1;
	  }
//-----------------------------------------------------------------------------
// calculate chi2 of the found seed
//-----------------------------------------------------------------------------
	  seed->fChi2All = 0;
	  for (int face=0; face<kNFaces; face++) {
	    int nh = seed->NHits(face);
	    for (int ih=0; ih<nh; ih++) {
	      const HitData_t* hd = seed->HitData(face,ih);

	      // const StrawHitPosition* shp  = hd->fPos;
	      CLHEP::Hep3Vector       dxyz = hd->fHit->posCLHEP()-seed->CofM; //shp->posCLHEP()-seed->CofM; // distance from hit to the center-of-gravity
//-----------------------------------------------------------------------------
// split into wire parallel and perpendicular components
//-----------------------------------------------------------------------------
	      const CLHEP::Hep3Vector& wdir = hd->fHit->wdirCLHEP();//fStraw->getDirection();
	      CLHEP::Hep3Vector d_par       = (dxyz.dot(wdir))/(wdir.dot(wdir))*wdir; 
	      CLHEP::Hep3Vector d_perp_z    = dxyz-d_par;
	      float  d_perp                 = d_perp_z.perp();
	      double sigw                   = hd->fSigW;
	      float  chi2_par               = (d_par.mag()/sigw)*(d_par.mag()/sigw);
	      float  chi2_perp              = (d_perp/_sigmaR)*(d_perp/_sigmaR);
	      float  chi2                   = chi2_par + chi2_perp;
	      seed->fChi2All               += chi2;
	    }
	  }
	  seed->fChi2All = seed->fChi2All/seed->fNHitsTot;
	}
//-----------------------------------------------------------------------------
// prune list of found seeds
//-----------------------------------------------------------------------------
	pruneSeeds(s);
      }
    }
  }

  // unflagging preseed hits if seed is not completed?
  // start with object, fill in as it goes, then decide whether or not to keep it

  // move to other faces, use phi of preseed pos and phi of panels to determine which panel to check
  // loop through hits with time constraint
  // for chi2, use distance along the wire/resolution and distance across the wire/TBD(1 cm?)
  // find a way to accommodate for preseed size

  // Find "average straw" of all candidates in both directions, then find intersection of these two average straws

//-----------------------------------------------------------------------------
// consider only good seeds
//-----------------------------------------------------------------------------
  void DeltaFinder::connectSeeds() {

    //loop over stations
    //start with a seed
    //move to next station
    //loop over seeds, look for one with similar xy and time (split into components impractical?)
    //if multiple, select one with best chi2
    //continue, incrementing over stations until all finished
    //what to do if there are gaps in the path?
    //update center of mass?

    for (int s=0; s<kNStations; ++s) {
      int pssize = _data.seedHolder[s].size();
      for (int ps=0; ps<pssize; ++ps) {
	DeltaSeed* seed = _data.seedHolder[s][ps];
//-----------------------------------------------------------------------------
// create new delta candidate if a seed has >= _minNFacesWithHits
//-----------------------------------------------------------------------------
	if (seed->fGood < 0)                                continue;
	if (seed->Used()   )                                continue;
	if (seed->fNFacesWithHits < _minNFacesWithHits)     continue;

	DeltaCandidate delta;
	delta.seed[s]       = seed;
	delta.CofM          = seed->CofM;
	delta.n_seeds       = 1;
	delta.fFirstStation = s;
	delta.fLastStation  = s;
	delta.fNHits        = seed->fNHitsTot;
	delta.fT0Min[s]     = seed->T0Min();
	delta.fT0Max[s]     = seed->T0Max();
	//	double t0           = (seed->fMaxTime+seed->fMinTime-_maxDriftTime)/2;
	//	delta.fTzSums.addPoint(seed->CofM.z(),t0);
//-----------------------------------------------------------------------------
// stations 6 and 13 are empty - account for that
//-----------------------------------------------------------------------------
	int sdist = 0;
	for (int s2=s+1; s2<kNStations; ++s2) {
	  if ((s2 == 6) || (s2 == 13)) sdist += 1;
	  int gap = s2-delta.fLastStation-1;
	  if ( gap > _maxGap+sdist) break;
//------------------------------------------------------------------------------
// never extend 1-seg candidates over 2 empty stations
//-----------------------------------------------------------------------------
	  if ((s2-delta.fLastStation == delta.fFirstStation) && (gap > 1)) break;
//-----------------------------------------------------------------------------
// predict T0
//-----------------------------------------------------------------------------
	  float t0max  = delta.fT0Max[delta.fLastStation] + _maxDtDs*(s2-delta.fLastStation);
	  float t0min  = delta.fT0Min[delta.fLastStation] - _maxDtDs*(s2-delta.fLastStation);
//-----------------------------------------------------------------------------
// find the closest seed in station s2
//-----------------------------------------------------------------------------
	  DeltaSeed* closest(NULL);
	  float      dxy, dxy_min(_maxDxy);

	  int ps2size = _data.seedHolder[s2].size();
	  for (int ps2=0; ps2<ps2size; ++ps2) {
	    DeltaSeed* seed2 = _data.seedHolder[s2][ps2];
	    if (seed2->fGood < 0)                            continue;
	    if (seed2->Used()   )                            continue;
	    if (seed2->fNFacesWithHits < _minNFacesWithHits) continue;
	    if (seed2->T0Max() - t0min < -10.)               continue; // *FIXME* make a parameter
	    if (seed2->T0Min() - t0max >  10.)               continue;
//-----------------------------------------------------------------------------
// seed2 T0 is consistent with the predicted T0
//-----------------------------------------------------------------------------
	    CLHEP::Hep3Vector dxyz = seed2->CofM-delta.CofM;
	    dxy                    = dxyz.perp();

	    if (dxy < dxy_min) {
	      closest = seed2;
	      dxy_min = dxy;
	    }
	  }

	  if (closest) {
	    delta.fLastStation = s2;
	    delta.dxy    [s2]  = dxy_min;
	    delta.seed   [s2]  = closest;
	    delta.CofM         = (delta.CofM*delta.n_seeds+closest->CofM)/(delta.n_seeds+1);
	    delta.n_seeds     += 1;
	    delta.fNHits      += closest->fNHitsTot;
	    sdist              = 0;
//-----------------------------------------------------------------------------
// redefine T0 limits (in the last station)
//-----------------------------------------------------------------------------
	    if (closest->T0Min() < t0min) delta.fT0Min[s2] = t0min;
	    else                          delta.fT0Min[s2] = closest->T0Min();

	    if (closest->T0Max() > t0max) delta.fT0Max[s2] = t0max;
	    else                          delta.fT0Max[s2] = closest->T0Max();
	  }
	}
//-----------------------------------------------------------------------------
// store only delta candidates with more than 2 stations
// for each station define expected T0min and T0max
//-----------------------------------------------------------------------------
	if (delta.n_seeds >= _minNSeeds) {
	  delta.fNumber = _data.deltaCandidateHolder.size();
	  int last = -1;
	  for (int station=delta.fFirstStation; station<=delta.fLastStation; station++) {
	    DeltaSeed* ds = delta.seed[station];
	    if (ds != NULL) { 
	      ds->fDeltaIndex = delta.fNumber;

	      if (last == -1) {
		delta.fT0Min[station] = ds->T0Min();
		delta.fT0Max[station] = ds->T0Max();
		last                  = station;
	      }
	      else {
		float t0min = delta.fT0Min[last]-_maxDtDs*(station-last);
		float t0max = delta.fT0Max[last]+_maxDtDs*(station-last);

		if (t0min < ds->T0Min()) delta.fT0Min[station] = ds->T0Min();
		else                     delta.fT0Min[station] = t0min;

		if (t0max > ds->T0Max()) delta.fT0Max[station] = ds->T0Min();
		else                     delta.fT0Max[station] = t0max;
	      }
	    }
	  }
	  delta.phi     = delta.CofM.phi();                  // calculate just once
	  _data.deltaCandidateHolder.push_back(delta);
	}
      }
    }
  }
}
//-----------------------------------------------------------------------------
// magic that makes this class a module.
//-----------------------------------------------------------------------------
DEFINE_ART_MODULE(mu2e::DeltaFinder)
//-----------------------------------------------------------------------------
// done
//-----------------------------------------------------------------------------
