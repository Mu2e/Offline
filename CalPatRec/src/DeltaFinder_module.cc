/////////////////////////////////////////////////////////////////////////////
// framework
//
// parameter defaults: CalPatRec/fcl/prolog.fcl
//////////////////////////////////////////////////////////////////////////////
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art_root_io/TFileService.h"

#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"

// conditions
#include "Offline/ConditionsService/inc/ConditionsHandle.hh"

// data
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitPosition.hh"
#include "Offline/RecoDataProducts/inc/StereoHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"

// diagnostics

#include "Offline/CalPatRec/inc/DeltaFinder_types.hh"

// #include "CalPatRec/inc/LsqSums2.hh"
#include "Offline/Mu2eUtilities/inc/ModuleHistToolBase.hh"
#include "art/Utilities/make_tool.h"

#include <algorithm>
#include <cmath>
#include "CLHEP/Vector/ThreeVector.h"
#include "Offline/Mu2eUtilities/inc/TwoLinePCA.hh"
#include "Offline/Mu2eUtilities/inc/polyAtan2.hh"

using namespace std;

using CLHEP::Hep3Vector;

namespace mu2e {

  using namespace DeltaFinderTypes;

  class DeltaFinder: public art::EDProducer {
  public:

    struct Config {
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag>     shCollTag        {Name("shCollTag"        ), Comment("SComboHit collection Name"   ) };
      fhicl::Atom<art::InputTag>     chCollTag        {Name("chCollTag"        ), Comment("ComboHit collection Name"    ) };
      fhicl::Atom<art::InputTag>     chfCollTag       {Name("chfCollTag"       ), Comment("StrawHitFlag collection Name") };
      fhicl::Atom<art::InputTag>     sdmcCollTag      {Name("sdmcCollTag"      ), Comment("StrawDigiMC collection Name" ) };
      fhicl::Atom<art::InputTag>     tpeakCollTag     {Name("tpeakCollTag"     ), Comment("Time peak collection Name"   ) };
      fhicl::Atom<int>               useTimePeaks     {Name("useTimePeaks"     ), Comment("to use time peaks set to 1"  ) };
      fhicl::Atom<int>               debugLevel       {Name("debugLevel"       ), Comment("debug level"                 ) };
      fhicl::Atom<int>               diagLevel        {Name("diagLevel"        ), Comment("diag level"                  ) };
      fhicl::Atom<int>               printErrors      {Name("printErrors"      ), Comment("print errors"                ) };
      fhicl::Atom<float>             minCaloDt        {Name("minCaloDt"        ), Comment("min Calo Dt"                 ) };
      fhicl::Atom<float>             maxCaloDt        {Name("maxCaloDt"        ), Comment("max Calo Dt"                 ) };
      fhicl::Atom<float>             meanPitchAngle   {Name("meanPitchAngle"   ), Comment("mean pitch angle"            ) };
      fhicl::Atom<float>             minHitTime       {Name("minHitTime"       ), Comment("min hit time"                ) };
      fhicl::Atom<int>               minNFacesWithHits{Name("minNFacesWithHits"), Comment("min N faces with hits"       ) };
      fhicl::Atom<int>               minNSeeds        {Name("minNSeeds"        ), Comment("min N seeds in a delta cand" ) };
      fhicl::Atom<int>               minDeltaNHits    {Name("minDeltaNHits"    ), Comment("min N combo  hits in a delta") };
      fhicl::Atom<float>             maxEleHitEnergy  {Name("maxEleHitEnergy"  ), Comment("max electron hit energy"     ) };
      fhicl::Atom<float>             minimumTime      {Name("minimumTime"      ), Comment("minimum time"                ) };
      fhicl::Atom<float>             maximumTime      {Name("maximumTime"      ), Comment("maximum time"                ) };
      fhicl::Atom<float>             maxHitSeedDt     {Name("maxHitSeedDt"     ), Comment("max DT(hit-seed)"            ) };
      fhicl::Atom<float>             maxChi2Seed      {Name("maxChi2Seed"      ), Comment("max seed chi2 (stereo)"      ) };
      fhicl::Atom<float>             maxChi2Radial    {Name("maxChi2Radial"    ), Comment("max chi2 (radial)"           ) };
      fhicl::Atom<float>             maxChi2All       {Name("maxChi2All"       ), Comment("max chi2 (all)"              ) };
      fhicl::Atom<float>             maxChi2SeedDelta {Name("maxChi2SeedDelta" ), Comment("max chi2 (seed-delta)"       ) };
      fhicl::Atom<float>             seedRes          {Name("seedRes"          ), Comment("stereo seed resolution"      ) };
      fhicl::Atom<float>             maxDxy           {Name("maxDxy"           ), Comment("max Dxy"                     ) };
      fhicl::Atom<int>               maxGap           {Name("maxGap"           ), Comment("max Gap"                     ) };
      fhicl::Atom<float>             sigmaR           {Name("sigmaR"           ), Comment("sigmaR"                      ) };
      fhicl::Atom<float>             maxDriftTime     {Name("maxDriftTime"     ), Comment("maxDriftTime"                ) };
      fhicl::Atom<float>             maxSeedDt        {Name("maxSeedDt"        ), Comment("maxSeedDt"                   ) };
      fhicl::Atom<float>             maxHitDt         {Name("maxHitDt"         ), Comment("maxHitDt"                    ) };
      fhicl::Atom<float>             maxStrawDt       {Name("maxStrawDt"       ), Comment("max straw Dt"                ) };
      fhicl::Atom<float>             maxDtDs          {Name("maxDtDs"          ), Comment("max Dt/Dstation"             ) };
      fhicl::Atom<float>             maxDtDc          {Name("maxDtDc"          ), Comment("max deltaT between deltas"   ) };
      fhicl::Atom<int>               writeStrawHits   {Name("writeStrawHits"   ), Comment("if 1, write SCH coll"        ) };
      fhicl::Atom<int>               filter           {Name("filter"           ), Comment("if 1, write only nonDelta CH") };
      fhicl::Atom<int>               testOrder        {Name("testOrder"        ), Comment("if 1, test order"            ) };
      fhicl::Atom<bool>              testHitMask      {Name("testHitMask"      ), Comment("if true, test hit mask"      ) };
      fhicl::Sequence<std::string>   goodHitMask      {Name("goodHitMask"      ), Comment("good hit mask"               ) };
      fhicl::Sequence<std::string>   bkgHitMask       {Name("bkgHitMask"       ), Comment("background hit mask"         ) };

      fhicl::Table<DeltaFinderTypes::Config> diagPlugin{Name("diagPlugin"      ), Comment("Diag plugin") };
    };

  protected:
//-----------------------------------------------------------------------------
// talk-to parameters: input collections and algorithm parameters
//-----------------------------------------------------------------------------
    art::InputTag   _shCollTag;
    art::InputTag   _chCollTag;
    art::InputTag   _chfCollTag;
    art::InputTag   _sdmcCollTag;
    art::InputTag   _tpeakCollTag;

    int             _useTimePeaks;
    float           _minCaloDt;
    float           _maxCaloDt;
    float           _meanPitchAngle;
    float           _minHitTime;           // min hit time
    int             _minNFacesWithHits;    // per station per seed
    int             _minNSeeds;            // min number of seeds in the delta electron cluster
    int             _minDeltaNHits;        // min number of hits of a delta candidate
    float           _maxEleHitEnergy;      //
    float           _minT;
    float           _maxT;
    float           _maxHitSeedDt;         //
    float           _maxChi2Seed;          //
    float           _maxChi2Neighbor;      //
    float           _maxChi2Radial;        //
    float           _maxChi2All;           // max chi2/N of a seed
    float           _maxChi2SeedDelta;     // max delta-seed chi2 for adding a seed to a delta
    float           _seedRes;              //
    float           _maxDxy;
    int             _maxGap;
    float           _sigmaR;
    float           _sigmaR2;              // _sigmaR^2
    float           _maxDriftTime;
    float           _maxSeedDt;            // +/- SeedDt is the time window for checking duplicate seeds
    float           _maxHitDt;
    float           _maxStrawDt;
    float           _maxDtDs;              // low-P electron travel time between two stations
    float           _maxDtDc;              // max deltaT between two delta candiates
    int             _writeStrawHits;
    int             _filter;

    int             _debugLevel;
    int             _diagLevel;
    int             _printErrors;
    int             _testOrder;
    bool            _testHitMask;
    StrawHitFlag    _goodHitMask;
    StrawHitFlag    _bkgHitMask;

    std::unique_ptr<ModuleHistToolBase> _hmanager;
//-----------------------------------------------------------------------------
// cache event/geometry objects
//-----------------------------------------------------------------------------
    const StrawHitCollection*    _shColl ;

    const Tracker*               _tracker;
    const DiskCalorimeter*       _calorimeter;

    float                        _tdbuff; // following Dave - time division buffer

    DeltaFinderTypes::Data_t     _data;              // all data used
    int                          _testOrderPrinted;

    int                          _nComboHits;
    int                          _nStrawHits;

    vector<const ComboHit*>      _v;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  public:
    explicit DeltaFinder(const art::EDProducer::Table<Config>& config);

  private:

    bool         findData            (const art::Event&  Evt);
    int          checkDuplicates     (int Station,
                                      int Face1, const HitData_t* Hit1,
                                      int Face2, const HitData_t* Hit2);

    void         completeSeed        (DeltaSeed* Seed);
    void         connectSeeds        ();
    void         findSeeds           (int Station, int Face);
    void         findSeeds           ();
    void         initTimeCluster     (DeltaCandidate* Delta, TimeCluster* Tc);
    int          mergeDeltaCandidates();
    int          orderHits           ();
    void         pruneSeeds          (int Station);
    int          recoverMissingHits  ();
    int          recoverSeed         (DeltaCandidate* Delta, int LastStation, int Station);
    int          recoverStation      (DeltaCandidate* Delta, int LastStation, int Station, int UseUsedHits, int RecoverSeeds);
    void         runDeltaFinder      ();
    double       seedDeltaChi2       (DeltaSeed* Seed, DeltaCandidate* Delta);
//-----------------------------------------------------------------------------
// overloaded methods of the module class
//-----------------------------------------------------------------------------
    void         beginJob() override;
    void         beginRun(art::Run& ARun) override;
    void         endJob  () override;
    void         produce(art::Event& e) override;
  };

//-----------------------------------------------------------------------------
  DeltaFinder::DeltaFinder(const art::EDProducer::Table<Config>& config):
    art::EDProducer{config},
    _shCollTag             (config().shCollTag()        ),
    _chCollTag             (config().chCollTag()        ),
    _chfCollTag            (config().chfCollTag()       ),
    _sdmcCollTag           (config().sdmcCollTag()      ),
    _tpeakCollTag          (config().tpeakCollTag()     ),
    _useTimePeaks          (config().useTimePeaks()     ),
    _minCaloDt             (config().minCaloDt()        ),
    _maxCaloDt             (config().maxCaloDt()        ),
    _meanPitchAngle        (config().meanPitchAngle()   ),
    _minHitTime            (config().minHitTime()       ),
    _minNFacesWithHits     (config().minNFacesWithHits()),
    _minNSeeds             (config().minNSeeds()        ),
    _minDeltaNHits         (config().minDeltaNHits()    ),
    _maxEleHitEnergy       (config().maxEleHitEnergy()  ),
    _minT                  (config().minimumTime()      ), // nsec
    _maxT                  (config().maximumTime()      ), // nsec
    _maxHitSeedDt          (config().maxHitSeedDt()     ), // nsec
    _maxChi2Seed           (config().maxChi2Seed()      ),
    _maxChi2Radial         (config().maxChi2Radial()    ),
    _maxChi2All            (config().maxChi2All()       ),
    _maxChi2SeedDelta      (config().maxChi2SeedDelta() ),
    _seedRes               (config().seedRes()          ),
    _maxDxy                (config().maxDxy()           ),
    _maxGap                (config().maxGap()           ),
    _sigmaR                (config().sigmaR()           ),
    _maxDriftTime          (config().maxDriftTime()     ),
    _maxSeedDt             (config().maxSeedDt()        ),
    _maxHitDt              (config().maxHitDt()         ),
    _maxStrawDt            (config().maxStrawDt()       ),
    _maxDtDs               (config().maxDtDs()          ),
    _maxDtDc               (config().maxDtDc()          ),
    _writeStrawHits        (config().writeStrawHits()   ),
    _filter                (config().filter()           ),
    _debugLevel            (config().debugLevel()       ),
    _diagLevel             (config().diagLevel()        ),
    _printErrors           (config().printErrors()      ),
    _testOrder             (config().testOrder()        ),
    _testHitMask           (config().testHitMask()      ),
    _goodHitMask           (config().goodHitMask()      ),
    _bkgHitMask            (config().bkgHitMask()       )
  {
    consumesMany<ComboHitCollection>(); // Necessary because fillStrawHitIndices calls getManyByType.

    produces<StrawHitFlagCollection>("ComboHits");
    if (_writeStrawHits == 1) produces<StrawHitFlagCollection>("StrawHits");
    if (_filter             ) produces<ComboHitCollection>();

                                        // this is a list of delta-electron candidates
    produces<TimeClusterCollection>();

    _testOrderPrinted = 0;
    _tdbuff           = 80.;                 // mm ... about less than 1 ns
    _sigmaR2          = _sigmaR*_sigmaR;

    if (_diagLevel != 0) _hmanager = art::make_tool  <ModuleHistToolBase>(config().diagPlugin.get_PSet());
    else                 _hmanager = std::make_unique<ModuleHistToolBase>();

    _data.chCollTag      = _chCollTag;
    _data.chfCollTag     = _chfCollTag;
    _data.sdmcCollTag    = _sdmcCollTag;
    _data.meanPitchAngle = _meanPitchAngle;

  }

  //-----------------------------------------------------------------------------
  void DeltaFinder::beginJob() {
    if (_diagLevel > 0) {
      art::ServiceHandle<art::TFileService> tfs;
      _hmanager->bookHistograms(tfs);
    }
  }

  //-----------------------------------------------------------------------------
  void DeltaFinder::endJob() {
  }

//-----------------------------------------------------------------------------
// create a Z-ordered representation of the tracker
//-----------------------------------------------------------------------------
  void DeltaFinder::beginRun(art::Run& aRun) {


    _data.InitGeometry();

    // mu2e::GeomHandle<mu2e::Tracker> tHandle;
    // _tracker      = tHandle.get();
    // _data.tracker = _tracker;

    // mu2e::GeomHandle<mu2e::DiskCalorimeter> ch;
    // _calorimeter = ch.get();

    // ChannelID cx, co;

    // int       nDisks    = _calorimeter->nDisk();
    // double    disk_z[2] = {0};//given in the tracker frame

    // for (int i=0; i<nDisks; ++i){
    //   Hep3Vector gpos = _calorimeter->disk(i).geomInfo().origin();
    //   Hep3Vector tpos = _calorimeter->geomUtil().mu2eToTracker(gpos);
    //   disk_z[i] = tpos.z();
    // }
//-----------------------------------------------------------------------------
// define station Z coordinates
//-----------------------------------------------------------------------------
    // for (unsigned ipl=0; ipl<_tracker->nPlanes(); ipl += 2) {
    //   const Plane* p1 = &_tracker->getPlane(ipl);
    //   const Plane* p2 = &_tracker->getPlane(ipl+1);
    //   _data.stationZ[ipl/2] = (p1->origin().z()+p2->origin().z())/2;
    // }

//     float     z_tracker_center(0.);
//     int       nPlanesPerStation(2);
//     double    station_z(0);

//     for (unsigned planeId=0; planeId<_tracker->nPlanes(); planeId++) {
//       const Plane* pln = &_tracker->getPlane(planeId);
//       int  ist = planeId/nPlanesPerStation;
//       int  ipl = planeId % nPlanesPerStation;
// //-----------------------------------------------------------------------------
// // calculate the time-of-flight between the station and each calorimeter disk
// // for a typical Conversion Electron
// //-----------------------------------------------------------------------------
//       if (ipl == 0) {
//         station_z = pln->origin().z();
//       }
//       else {
//         station_z = (station_z + pln->origin().z())/2.;
//         for (int iDisk=0; iDisk<nDisks; ++iDisk){
//           _stationToCaloTOF[iDisk][ist] = (disk_z[iDisk] - station_z)/sin(_meanPitchAngle)/CLHEP::c_light;
//         }
//       }

//       for (unsigned ipn=0; ipn<pln->nPanels(); ipn++) {
//         const Panel* panel = &pln->getPanel(ipn);
//         int face;
//         if (panel->id().getPanel() % 2 == 0) face = 0;
//         else                                 face = 1;
//         for (unsigned il=0; il<panel->nLayers(); ++il) {
//           cx.Station   = ist;
//           cx.Plane     = ipl;
//           cx.Face      = face;
//           cx.Panel     = ipn;
//           cx.Layer     = il;
//           orderID (&cx, &co);
//           int os       = co.Station;
//           int of       = co.Face;
//           int op       = co.Panel;
//           PanelZ_t* pz = &_data.oTracker[os][of][op];
//           pz->fPanel   = panel;
// //-----------------------------------------------------------------------------
// // panel caches phi of its center and the z, phi runs from 0 to 2*pi
// //-----------------------------------------------------------------------------
//           pz->wx  = panel->straw0Direction().x();
//           pz->wy  = panel->straw0Direction().y();
//           pz->phi = panel->straw0MidPoint().phi();
//           pz->z   = (panel->getStraw(0).getMidPoint().z()+panel->getStraw(1).getMidPoint().z())/2.;
//           int  uniqueFaceId = ipl*mu2e::StrawId::_nfaces + of;
//           _faceTOF[uniqueFaceId] = (z_tracker_center - pz->z)/sin(_meanPitchAngle)/CLHEP::c_light;
//         }
//       }
//       _data.stationUsed[ist] = 1;
//     }
//-----------------------------------------------------------------------------
// it is enough to print that once
//-----------------------------------------------------------------------------
    if (_testOrder && (_testOrderPrinted == 0)) {
      _data.testOrderID  ();
      _data.testdeOrderID();
      _testOrderPrinted = 1;
    }

    if (_diagLevel != 0) _hmanager->debug(&_data,1);
  }

//-----------------------------------------------------------------------------
// make sure the two hits used to make a new seed are not a part of an already found seed
//-----------------------------------------------------------------------------
  int DeltaFinder::checkDuplicates(int Station, int Face1, const HitData_t* Hit1, int Face2, const HitData_t* Hit2) {

    int rc(0);

    int nseeds = _data.NSeeds(Station);
    for (int i=0; i<nseeds; i++) {
      DeltaSeed* seed = _data.deltaSeed(Station,i);

      if ((seed->hitData[Face1] == Hit1) and (seed->hitData[Face2] == Hit2)) {
//-----------------------------------------------------------------------------
// 'seed' contains both Hit1 and Hit2, done
//-----------------------------------------------------------------------------
        rc = 1;
                                                                          break;
      }
    }
    return rc;
  }

//-----------------------------------------------------------------------------
// input 'Seed' has two hits , non parallel wires
// try to add more close hits to it (one hit per face)
//-----------------------------------------------------------------------------
  void DeltaFinder::completeSeed(DeltaSeed* Seed) {

    assert ((Seed->SFace(1) >= 0) and (Seed->NHits() == 2));
//-----------------------------------------------------------------------------
// loop over remaining faces, 'f2' - face in question
//-----------------------------------------------------------------------------
    int station = Seed->Station();

    float xseed = Seed->CofM.x();
    float yseed = Seed->CofM.y();

    for (int face=0; face<kNFaces; face++) {
      if (Seed->fFaceProcessed[face] == 1)                            continue;
//-----------------------------------------------------------------------------
// face is different from the two first faces used
//-----------------------------------------------------------------------------
      for (int p2=0; p2<3; ++p2) {
        PanelZ_t* panelz = &_data.oTracker[station][face][p2];
        double dphi      = Seed->Phi()-panelz->phi;
        if (dphi < -M_PI) dphi += 2*M_PI;
        if (dphi >  M_PI) dphi -= 2*M_PI;
        if (fabs(dphi) >= M_PI/3)                                     continue;
//-----------------------------------------------------------------------------
// panel overlaps with the seed, look at its hits
//-----------------------------------------------------------------------------
        double nx    = panelz->wx;
        double ny    = panelz->wy;

        float sxy_dot_w = xseed*nx+yseed*ny;

        assert(sxy_dot_w != 1.e10);             // test
        if (sxy_dot_w > 1.e5) {
          printf("emoe\n"); // test
        }

        int    nhits = panelz->fHitData->size();
        for (int ih=0; ih<nhits; ih++) {
//-----------------------------------------------------------------------------
// 2017-10-05 PM: consider all hits
// hit time should be consistent with the already existing times - the difference
// between any two measured hit times should not exceed _maxDriftTime
// (_maxDriftTime represents the maximal drift time in the straw, should there be some tolerance?)
//-----------------------------------------------------------------------------
          HitData_t* hd            = &(*panelz->fHitData)[ih];
          float corr_time          = hd->fCorrTime;

          if (corr_time-Seed->T0Max() > _maxHitSeedDt)                break;
          if (Seed->T0Min()-corr_time > _maxHitSeedDt)                continue;

          const ComboHit*   ch     = hd->fHit;
          double dx                = ch->pos().x()-xseed;
          double dy                = ch->pos().y()-yseed;
//-----------------------------------------------------------------------------
// split into wire parallel and perpendicular components
// assume wire direction vector is normalized to 1
//-----------------------------------------------------------------------------
          float dxy_dot_wdir       = dx*nx+dy*ny;

          float dx_perp            = dx-dxy_dot_wdir*nx;
          float dy_perp            = dy-dxy_dot_wdir*ny;
          float d_perp_2           = dx_perp*dx_perp+dy_perp*dy_perp;

          // add an uncertainty on the intersection, can do better :

          float chi2_par           = (dxy_dot_wdir*dxy_dot_wdir)/(hd->fSigW2+_seedRes*_seedRes);
          float chi2_perp          = d_perp_2/_sigmaR2;
          float chi2               = chi2_par + chi2_perp;
          if (chi2 >= _maxChi2Radial)                                 continue;
//-----------------------------------------------------------------------------
// add hit
//-----------------------------------------------------------------------------
          hd->fChi2Min             = chi2;
          Seed->AddHit(hd,face);
        }
      }
    }
//-----------------------------------------------------------------------------
// update seed time and X and Y coordinates, accurate knowledge of Z is not very relevant
//-----------------------------------------------------------------------------
    Seed->CalculateCogAndChi2(_sigmaR2);
  }

//-----------------------------------------------------------------------------
// loop over stations
// start with a seed
// move to next station
// loop over seeds, look for one with similar xy and time (split into components impractical?)
// if multiple, select one with best chi2
// continue, incrementing over stations until all finished
// what to do if there are gaps in the path?
// update center of mass?
//-----------------------------------------------------------------------------
  void DeltaFinder::connectSeeds() {

    for (int is=0; is<kNStations; is++) {
//-----------------------------------------------------------------------------
// 1. loop over existing seeds and match them to existing delta candidates
//-----------------------------------------------------------------------------
      int nseeds = _data.NSeeds(is);
      for (int ids=0; ids<nseeds; ids++) {
        DeltaSeed* seed = _data.deltaSeed(is,ids);

        if (! seed->Good() )                                          continue;
        if (seed->Used()   )                                          continue;
        //        if (seed->fNFacesWithHits < _minNFacesWithHits)               continue;
//-----------------------------------------------------------------------------
// first, loop over existing delta candidates and try to associate the seed
// with one of them
//-----------------------------------------------------------------------------
        DeltaCandidate* closest(nullptr);
        float           chi2min (_maxChi2SeedDelta);

        int ndelta = _data.listOfDeltaCandidates.size();
        for (int idc=0; idc<ndelta; idc++) {
          DeltaCandidate* dc = _data.deltaCandidate(idc);
//-----------------------------------------------------------------------------
// skip candidates already merged with others
//-----------------------------------------------------------------------------
          if (dc->Active() == 0      )                                continue;
          if (dc->seed[is] != nullptr)                                continue;
          int last = dc->LastStation();
//-----------------------------------------------------------------------------
// a delta candidate starting from a seed in a previous station may already have
// a seed in this station found, skip such a candidate
//-----------------------------------------------------------------------------
          if (last == is)                                             continue;
          int gap  = is-last;
          if (gap > _maxGap)                                          continue;
          float t0 = dc->T0(is);
          assert(t0 > 0);
          float dt = t0-seed->TMean();
          if (fabs(dt) > _maxSeedDt+_maxDtDs*gap)                     continue;
//-----------------------------------------------------------------------------
// the time is OK'ish - checks should be implemented more accurately (FIXME)
// look at the coordinates
//-----------------------------------------------------------------------------
          float chi2 = seedDeltaChi2(seed,dc);
          if (chi2 < chi2min) {
//-----------------------------------------------------------------------------
// everything matches, new closest delta candidate
//-----------------------------------------------------------------------------
            closest = dc;
            chi2min = chi2;
          }
        }
//-----------------------------------------------------------------------------
// if a DeltaSeed has been "attached" to a DeltaCandidate, this is it.
//-----------------------------------------------------------------------------
        if (closest) {
          closest->AddSeed(seed,is);
          seed->fChi2Delta = chi2min;
                                                                      continue;
        }
//-----------------------------------------------------------------------------
// DeltaSeed has not been linked to any existing delta candidate, create
// a new delta candidate and see if it is good enough
//-----------------------------------------------------------------------------
        DeltaCandidate delta(ndelta,seed,is);
//-----------------------------------------------------------------------------
// first, try to extend it backwards, in a compact way, to pick missing single hits
// seeds should've already been picked up !
//-----------------------------------------------------------------------------
        for (int is2=is-1; is2>=0; is2--) {
          recoverStation(&delta,delta.fFirstStation,is2,1,0);
//-----------------------------------------------------------------------------
// continue only if found something, allow one gap
//-----------------------------------------------------------------------------
          if (delta.fFirstStation-is2 > 2) break;
        }
//-----------------------------------------------------------------------------
// next, try to extend it forward, by one step, use of seeds is allowed
//-----------------------------------------------------------------------------
        if (is < kNStations-1) {
          recoverStation(&delta,is,is+1,1,1);
        }
//-----------------------------------------------------------------------------
// store only delta candidates with hits in more than 2 stations
// for each station define expected T0min and T0max
// to keep _minNSeeds=2 need to look forward...
//-----------------------------------------------------------------------------
        if (delta.fNSeeds >= _minNSeeds) {
          _data.listOfDeltaCandidates.push_back(delta);
        }
        else {
//-----------------------------------------------------------------------------
// mark all seeds as unassigned, this should be happening only in 1-seed case
//-----------------------------------------------------------------------------
          for (int is=delta.fFirstStation; is<=delta.fLastStation; is++) {
            DeltaSeed* ds = delta.seed[is];
            if (ds) ds->fDeltaIndex = -1;
          }
        }
      }
//-----------------------------------------------------------------------------
// all seeds in a given station processed
// loop over existing delta candidates which do not have seeds in a given station
// and see if can pick up single hits
//-----------------------------------------------------------------------------
      int ndelta = _data.listOfDeltaCandidates.size();
      for (int idc=0; idc<ndelta; idc++) {
        DeltaCandidate* dc = _data.deltaCandidate(idc);
        int last = dc->LastStation();
        if (last != is-1)                                           continue;
//-----------------------------------------------------------------------------
// if a delta candidate has been created in this routine, time limits
// may not be defined. Make sure they are
//-----------------------------------------------------------------------------
        recoverStation(dc,last,is,1,0);
      }
    }

//-----------------------------------------------------------------------------
// at this point we have a set of delta candidates, which might need to be merged
//-----------------------------------------------------------------------------
    mergeDeltaCandidates();
  }

//-----------------------------------------------------------------------------
  bool DeltaFinder::findData(const art::Event& Evt) {
    _data.chcol    = nullptr;
    _data.tpeakcol = nullptr;

    if (_useTimePeaks == 1){
      auto tpeakH = Evt.getValidHandle<mu2e::TimeClusterCollection>(_tpeakCollTag);
      _data.tpeakcol = tpeakH.product();
    }

    auto chcH   = Evt.getValidHandle<mu2e::ComboHitCollection>(_chCollTag);
    _data.chcol = chcH.product();

    auto chcfH    = Evt.getValidHandle<mu2e::StrawHitFlagCollection>(_chfCollTag);
    _data.chfColl = chcfH.product();

    auto shcH = Evt.getValidHandle<mu2e::StrawHitCollection>(_shCollTag);
    _shColl   = shcH.product();

    return (_data.chcol != nullptr) and (_data.chfColl != nullptr) and (_shColl != nullptr);
  }

//-----------------------------------------------------------------------------
// find delta electron seeds in 'Station' with hits in faces 'f' and 'f+1'
// do not consider proton hits with eDep > _minHitEnergy
//-----------------------------------------------------------------------------
  void DeltaFinder::findSeeds(int Station, int Face) {

    for (int p=0; p<3; ++p) {                        // loop over panels in this face
      PanelZ_t* panelz1 = &_data.oTracker[Station][Face][p];
      double nx1        = panelz1->wx;
      double ny1        = panelz1->wy;
      int    nh1        = panelz1->fHitData->size();
      for (int h1=0; h1<nh1; ++h1) {
//-----------------------------------------------------------------------------
// hit has not been used yet to start a seed, however it could've been used as a second seed
//-----------------------------------------------------------------------------
        HitData_t*      hd1 = &(*panelz1->fHitData)[h1];
        float           ct1 = hd1->fCorrTime;
        double x1           = hd1->fHit->pos().x();
        double y1           = hd1->fHit->pos().y();

        int counter         = 0;                // number of stereo candidate hits close to set up counter
//-----------------------------------------------------------------------------
// loop over the second faces
//-----------------------------------------------------------------------------
        for (int f2=Face+1; f2<kNFaces; f2++) {
          for (int p2=0; p2<3; ++p2) {         // loop over panels
            PanelZ_t* panelz2 = &_data.oTracker[Station][f2][p2];
//-----------------------------------------------------------------------------
// check if the two panels overlap in XY
// 2D angle between the vectors pointing to the panel centers, can't be greater than 2*pi/3
//-----------------------------------------------------------------------------
            float dphi = panelz2->phi - panelz1->phi;
            if (dphi < -M_PI) dphi += 2*M_PI;
            if (dphi >  M_PI) dphi -= 2*M_PI;
            if (abs(dphi) >= 2*M_PI/3.)                                continue;
//-----------------------------------------------------------------------------
// panels do overlap, check the time. tmin and tmax also detect panels w/o hits
//-----------------------------------------------------------------------------
            if (panelz2->tmin - ct1 > _maxDriftTime)                   continue;
            if (ct1 - panelz2->tmax > _maxDriftTime)                   continue;

            double nx2   = panelz2->wx;
            double ny2   = panelz2->wy;
            double n1n2  = nx1*nx2+ny1*ny2;
            double q12   = 1-n1n2*n1n2;
            double res_z = (panelz1->z+panelz2->z)/2;

            int    nh2 = panelz2->fHitData->size();
            for (int h2=0; h2<nh2;++h2) {
              HitData_t* hd2 = &(*panelz2->fHitData)[h2];
              float      ct2 = hd2->fCorrTime;
//-----------------------------------------------------------------------------
// hits are ordered in time, so if ct2-ct > _maxDriftTime, can proceed with the next panel
//-----------------------------------------------------------------------------
              if (ct2 - ct1 > _maxDriftTime)                           break;
              if (ct1 - ct2 > _maxDriftTime)                           continue;
              ++counter;                                            // number of hits close to the first one
//-----------------------------------------------------------------------------
// intersect the two straws, we need coordinates of the intersection point and
// two distances from hits to the intersection point, 4 numbers in total
//-----------------------------------------------------------------------------
              double x2    = hd2->fHit->pos().x();
              double y2    = hd2->fHit->pos().y();

              double r12n1 = (x1-x2)*nx1+(y1-y2)*ny1;
              double r12n2 = (x1-x2)*nx2+(y1-y2)*ny2;

              double wd1   = -(r12n2*n1n2-r12n1)/q12;

              double res_x = x1-nx1*wd1;
              double res_y = y1-ny1*wd1;

              double wd2   = -(r12n2-n1n2*r12n1)/q12;
//-----------------------------------------------------------------------------
// require both hits to be close enough to the intersection point
//-----------------------------------------------------------------------------
              float hd1_chi2 = wd1*wd1/hd1->fSigW2;
              float hd2_chi2 = wd2*wd2/hd2->fSigW2;
              if (hd1_chi2 > _maxChi2Seed)                            continue;
              if (hd2_chi2 > _maxChi2Seed)                            continue;
              if ((hd1_chi2+hd2_chi2)/2 > _maxChi2Seed)               continue;
//-----------------------------------------------------------------------------
// check whether there already is a seed containing both hits
//-----------------------------------------------------------------------------
              int is_duplicate  = checkDuplicates(Station,Face,hd1,f2,hd2);
              if (is_duplicate)                               continue;
//-----------------------------------------------------------------------------
// new seed : an intersection of two wires coresponsing to close in time combo hits
//-----------------------------------------------------------------------------
              hd1->fChi2Min     = hd1_chi2;
              hd2->fChi2Min     = hd2_chi2;

              DeltaSeed* seed   = _data.NewDeltaSeed(Station,Face,hd1,f2,hd2);

              seed->CofM.SetXYZ(res_x,res_y,res_z);
              seed->fPhi = polyAtan2(res_y,res_x);
//-----------------------------------------------------------------------------
// mark both hits as a part of a seed, so they would not be used individually
// - see HitData_t::Used()
//-----------------------------------------------------------------------------
              hd1->fSeed  = seed;
              hd2->fSeed  = seed;
//-----------------------------------------------------------------------------
// complete search for hits of this seed, mark it BAD if a proton
// in principle, could place "high-charge" seeds into a separate list
// that should improve the performance
// make eDep cut value a module parameter
//-----------------------------------------------------------------------------
              completeSeed(seed);
              if (seed->EDep() > 0.005) {
                seed->fGood = -2000-seed->fIndex;
              }
            }
          }
        }
 //-----------------------------------------------------------------------------
 // this is needed for diagnostics only
 //-----------------------------------------------------------------------------
        // if (_diagLevel > 0) {
        //   hd1->fNSecondHits  = counter ;
        // }
      }
    }
  }

//-----------------------------------------------------------------------------
// TODO: update the time as more hits are added
//-----------------------------------------------------------------------------
  void DeltaFinder::findSeeds() {

    for (int s=0; s<kNStations; ++s) {
      for (int face=0; face<kNFaces-1; face++) {
        findSeeds(s,face);
      }
      pruneSeeds(s);
    }
  }

//-----------------------------------------------------------------------------
// define the time cluster parameters starting from a DeltaCandidate
//-----------------------------------------------------------------------------
  void DeltaFinder::initTimeCluster(DeltaCandidate* Dc, TimeCluster* Tc) {
  }

//-----------------------------------------------------------------------------
// merge Delta Candidates : check for duplicates !
//-----------------------------------------------------------------------------
  int DeltaFinder::mergeDeltaCandidates() {
    int rc(0);
    float max_d2(20*20);  // mm^2, to be adjusted FIXME

    int ndelta = _data.listOfDeltaCandidates.size();

    for (int i1=0; i1<ndelta-1; i1++) {
      DeltaCandidate* dc1 = _data.deltaCandidate(i1);
      if (dc1->Active() == 0)                                          continue;
      float x1 = dc1->CofM.x();
      float y1 = dc1->CofM.y();
      for (int i2=i1+1; i2<ndelta; i2++) {
        DeltaCandidate* dc2 = _data.deltaCandidate(i2);
        if (dc2->Active() == 0)                                        continue;
        float x2 = dc2->CofM.x();
        float y2 = dc2->CofM.y();
        float d2 = (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);
        if (d2 > max_d2)                                               continue;
//-----------------------------------------------------------------------------
// too lazy to extrapolate the time to the same Z...  ** FIXME
//-----------------------------------------------------------------------------
        float t1 = dc1->T0(dc1->LastStation());
        float t2 = dc2->T0(dc2->FirstStation());
//-----------------------------------------------------------------------------
// time check could be done more intelligently - compare time at the same Z
//-----------------------------------------------------------------------------
        if (fabs(t1-t2) > _maxDtDc)                                    continue;

        int dds = dc2->FirstStation()-dc1->LastStation();
        if (dds < 0) {
          if (_printErrors) {
            printf("ERROR in DeltaFinder::%s:",__func__);
            printf("i1, i2, dc1->LastStation, dc2->FirstStation: %2i %2i %2i %2i \n",
                   i1,i2,dc1->LastStation(),dc2->FirstStation());
          }
                                                                       continue;
        }
        else if (dds > _maxGap)                                        continue;
//-----------------------------------------------------------------------------
// merge two delta candidates, not too far from each other in Z,
// leave dc1 active, mark dc2 as not active
//-----------------------------------------------------------------------------
        dc1->MergeDeltaCandidate(dc2,_printErrors);
        dc2->SetIndex(-1000-dc1->Index());
      }
    }
    return rc;
  }

//-----------------------------------------------------------------------------
// Custom comparator to sort in ascending order
//-----------------------------------------------------------------------------
  bool comparator(const ComboHit*& a, const ComboHit*& b) {
    return a->correctedTime() < b->correctedTime();
  }

//------------------------------------------------------------------------------
// I'd love to use the hit flags, however that is confusing:
// - hits with very large deltaT get placed to the middle of the wire and not flagged,
// - however, some hits within the fiducial get flagged with the ::radsel flag...
// use only "good" hits
//-----------------------------------------------------------------------------

  int DeltaFinder::orderHits() {
    ChannelID cx, co;

//-----------------------------------------------------------------------------
// vector of pointers to thecombo hits, ordered in time. Initial list is not touched
//-----------------------------------------------------------------------------
    _v.clear();
    _v.resize(_nComboHits);

    for (int i=0; i<_nComboHits; i++) {
      // _v.push_back(&(*_data.chcol)[i]); //  = &(*_chcol)[i];
      _v[i] = &(*_data.chcol)[i];
    }

    std::sort(_v.begin(), _v.end(), comparator);
//-----------------------------------------------------------------------------
// at this point all hits are ordered in time
//-----------------------------------------------------------------------------
    const ComboHit* ch0 = &(*_data.chcol)[0];

    for (int ih=0; ih<_nComboHits; ih++) {
      const ComboHit* ch = _v[ih];

      int loc = ch-ch0;
      const StrawHitFlag* flag   = &(*_data.chfColl)[loc];
      if (_testHitMask && (! flag->hasAllProperties(_goodHitMask) || flag->hasAnyProperty(_bkgHitMask)) ) continue;

      float corr_time    = ch->correctedTime();

      // if (ch->energyDep() > _maxEleHitEnergy             )  continue;
      if ((corr_time      < _minT) || (corr_time > _maxT))  continue;

      cx.Station                 = ch->strawId().station();
      cx.Plane                   = ch->strawId().plane() % 2;
      cx.Face                    = -1;
      cx.Panel                   = ch->strawId().panel();
      cx.Layer                   = ch->strawId().layer();

                                              // get Z-ordered location
      Data_t::orderID(&cx, &co);

      int os       = co.Station;
      int of       = co.Face;
      int op       = co.Panel;
      // int ol       = co.Layer;

      if (_useTimePeaks == 1) {
        bool               intime(false);
        int                nTPeaks  = _data.tpeakcol->size();
        const CaloCluster* cl(nullptr);
        int                iDisk(-1);

        for (int i=0; i<nTPeaks; ++i) {
          cl    = _data.tpeakcol->at(i).caloCluster().get();
          if (cl == nullptr) {
            printf(">>> DeltaFinder::orderHits() no CaloCluster found within the time peak %i\n", i);
            continue;
          }
          iDisk = cl->diskID();
          double    dt = cl->time() - (corr_time + _data.stationToCaloTOF[iDisk][os]);
          if ( (dt < _maxCaloDt) && (dt > _minCaloDt) ) {
            intime = true;
            break;
          }
        }
        if (!intime)                                    continue;
      }

      PanelZ_t* pz = &_data.oTracker[os][of][op];

      if (_printErrors) {
        if ((os < 0) || (os >= kNStations     )) printf(" >>> ERROR: wrong station number: %i\n",os);
        if ((of < 0) || (of >= kNFaces        )) printf(" >>> ERROR: wrong face    number: %i\n",of);
        if ((op < 0) || (op >= kNPanelsPerFace)) printf(" >>> ERROR: wrong panel   number: %i\n",op);
        // if ((ol < 0) || (ol >= 2              )) printf(" >>> ERROR: wrong layer   number: %i\n",ol);
      }

      pz->fHitData->push_back(HitData_t(ch));

      if (pz->tmin > corr_time) pz->tmin = corr_time;
      if (pz->tmax < corr_time) pz->tmax = corr_time;
    }

    return 0;
  }

//-----------------------------------------------------------------------------
  void DeltaFinder::produce(art::Event& Event) {

    if (_debugLevel) printf(">>> DeltaFinder::produce  event number: %10i\n",Event.event());
//-----------------------------------------------------------------------------
// clear memory in the beginning of event processing and cache event pointer
//-----------------------------------------------------------------------------
    _data.InitEvent(&Event,_debugLevel);
//-----------------------------------------------------------------------------
// process event
//-----------------------------------------------------------------------------
    if (! findData(Event)) {
      const char* message = "mu2e::DeltaFinder_module::produce: data missing or incomplete";
      throw cet::exception("RECO")<< message << endl;
    }

    _nComboHits = _data.chcol->size();
    _nStrawHits = _shColl->size();

    runDeltaFinder();
//-----------------------------------------------------------------------------
// form output - flag combo hits
//-----------------------------------------------------------------------------
    unique_ptr<StrawHitFlagCollection> up_chfcol(new StrawHitFlagCollection(_nComboHits));
    _data.outputChfColl = up_chfcol.get();

    for (int i=0; i<_nComboHits; i++) {
      StrawHitFlag flag;
      const ComboHit* ch = &(*_data.chcol)[i];
      int ind            = ch->indexArray().at(0);
      const StrawHit* sh = &_shColl->at(ind);
//-----------------------------------------------------------------------------
// make decision based on the first straw hit
//-----------------------------------------------------------------------------
      if (sh->energyDep() < _maxEleHitEnergy) flag.merge(StrawHitFlag::energysel);
      (*_data.outputChfColl)[i] = flag;
    }

    const ComboHit* ch0(0);
    if (_nComboHits > 0) ch0 = &_data.chcol->at(0);

    StrawHitFlag deltamask(StrawHitFlag::bkg);

    unique_ptr<TimeClusterCollection>  tcColl(new TimeClusterCollection);

    int ndeltas = _data.listOfDeltaCandidates.size();

    for (int i=0; i<ndeltas; i++) {
      DeltaCandidate* dc = &_data.listOfDeltaCandidates.at(i);
//-----------------------------------------------------------------------------
// skip merged in delta candidates
// also require a delta candidate to have at least 5 hits
// do not consider proton stub candidates (those with <EDep> > 0.004)
//-----------------------------------------------------------------------------
      if (dc->Active() == 0)                                          continue;
      if (dc->NHits () < _minDeltaNHits)                              continue;
      if (dc->EDep  () > 0.004         )                              continue;
      for (int station=dc->fFirstStation; station<=dc->fLastStation; station++) {
        DeltaSeed* ds = dc->seed[station];
        if (ds != nullptr) {
//-----------------------------------------------------------------------------
// loop over the hits and flag each of them as delta
//-----------------------------------------------------------------------------
          for (int face=0; face<kNFaces; face++) {
            const HitData_t* hd = ds->HitData(face);
            if (hd == nullptr)                                        continue;
            int loc = hd->fHit-ch0;
            _data.outputChfColl->at(loc).merge(deltamask);
          }
        }
      }
//-----------------------------------------------------------------------------
// make a time cluster out of each active DeltaCandidate
//-----------------------------------------------------------------------------
      TimeCluster new_tc;
      initTimeCluster(dc,&new_tc);
      tcColl->push_back(new_tc);
    }

    Event.put(std::move(tcColl));
//-----------------------------------------------------------------------------
// in the end of event processing fill diagnostic histograms
//-----------------------------------------------------------------------------
    if (_diagLevel  > 0) _hmanager->fillHistograms(&_data);
    if (_debugLevel > 0) _hmanager->debug(&_data,2);
//-----------------------------------------------------------------------------
// create the collection of StrawHitFlag for the StrawHitCollection
//-----------------------------------------------------------------------------
    if (_writeStrawHits == 1) {
                                        // first, copy over the original flags

      std::unique_ptr<StrawHitFlagCollection> shfcol(new StrawHitFlagCollection(_nStrawHits));

      for(int ich=0; ich<_nComboHits; ich++) {
        const ComboHit* ch   = &(*_data.chcol )[ich];
        StrawHitFlag    flag =  (*_data.outputChfColl)[ich];
        flag.merge(ch->flag());
        for (auto ish : ch->indexArray()) {
          (*shfcol)[ish] = flag;
        }
      }

      Event.put(std::move(shfcol),"StrawHits");
    }
//-----------------------------------------------------------------------------
// moving in the end, after diagnostics plugin routines have been called - move
// invalidates the original pointer...
//-----------------------------------------------------------------------------
    Event.put(std::move(up_chfcol),"ComboHits");
  }

//-----------------------------------------------------------------------------
// some of found seeds could be duplicates or ghosts
// in case two DeltaSeeds share the first seed hit, leave only the best one
// the seeds we're loooping over have been reconstructed within the same station
// also reject seeds with Chi2Tot > _maxChi2Tot=10
//-----------------------------------------------------------------------------
  void DeltaFinder::pruneSeeds(int Station) {

    int nseeds =  _data.NSeeds(Station);

    for (int i1=0; i1<nseeds-1; i1++) {
      DeltaSeed* ds1 = _data.deltaSeed(Station,i1);
      if (ds1->fGood < 0)                                             continue;

      if (ds1->Chi2AllN() > _maxChi2All) {
        ds1->fGood = -1000-i1;
                                                                      continue;
      }

      float tmean1 = ds1->TMean();

      for (int i2=i1+1; i2<nseeds; i2++) {
        DeltaSeed* ds2 = _data.deltaSeed(Station,i2);
        if (ds2->fGood < 0)                                           continue;

        if (ds2->Chi2AllN() > _maxChi2All) {
          ds2->fGood = -1000-i2;
                                                                      continue;
        }

        float tmean2 = ds2->TMean();

        if (fabs(tmean1-tmean2) > _maxSeedDt)                         continue;
//-----------------------------------------------------------------------------
// the two segments are close in time , both have acceptable chi2's
// *FIXME* didn't check distance !!!!!
// so far, allow duplicates during the search
// the two DeltaSeeds share could have significantly overlapping hit content
//-----------------------------------------------------------------------------
        int noverlap            = 0;
        int nfaces_with_overlap = 0;
        for (int face=0; face<kNFaces; face++) {
          const HitData_t* hh1 = ds1->hitData[face];
          const HitData_t* hh2 = ds2->hitData[face];
          if (hh1 and (hh1 == hh2)) {
            noverlap            += 1;
            nfaces_with_overlap += 1;
          }
        }

        if (nfaces_with_overlap > 1) {
//-----------------------------------------------------------------------------
// overlap significant, leave in only one DeltaSeed - which one?
//-----------------------------------------------------------------------------
          if (ds1->fNFacesWithHits > ds2->fNFacesWithHits) {
            ds2->fGood = -1000-i1;
          }
          else if (ds2->fNFacesWithHits > ds1->fNFacesWithHits) {
            ds1->fGood = -1000-i2;
            break;
          }
          else {
//-----------------------------------------------------------------------------
//both seeds have the same number of hits - compare chi2's
//-----------------------------------------------------------------------------
            if (ds1->Chi2AllN() <  ds2->Chi2AllN()) {
              ds2->fGood = -1000-i1;
            }
            else {
              ds1->fGood = -1000-i2;
              break;
            }
          }
        }
        else if (nfaces_with_overlap > 0) {
//-----------------------------------------------------------------------------
// only one overlapping hit
// special treatment of 2-hit seeds to reduce the number of ghosts
//-----------------------------------------------------------------------------
          if (ds1->fNFacesWithHits == 2) {
            if (ds2->fNFacesWithHits > 2) {
              ds1->fGood = -1000-i2;
              break;
            }
            else {
//-----------------------------------------------------------------------------
// the second seed also has 2 faces with hits
//-----------------------------------------------------------------------------
              if (ds1->Chi2AllN() <  ds2->Chi2AllN()) ds2->fGood = -1000-i1;
              else {
                ds1->fGood = -1000-i2;
                break;
              }
            }
          }
          else {
//-----------------------------------------------------------------------------
// the first seed has N>2 hits
//-----------------------------------------------------------------------------
            if (ds2->fNFacesWithHits == 2) {
//-----------------------------------------------------------------------------
// the 2nd seed has only 2 hits and there is an overlap
//-----------------------------------------------------------------------------
              ds2->fGood = -1000-i1;
            }
            else {
//-----------------------------------------------------------------------------
// the second seed also has N>2 faces with hits, but there is only one overlap
// leave both seeds in
//-----------------------------------------------------------------------------
            }
          }
        }
      }
    }
  }

//------------------------------------------------------------------------------
// start from looking at the "holes" in the seed pattern
// delta candidates in the list are already required to have at least 2 segments
// extend them outwards by one station
//-----------------------------------------------------------------------------
  int DeltaFinder::recoverMissingHits() {

    int ndelta = _data.listOfDeltaCandidates.size();
    for (int idelta=0; idelta<ndelta; idelta++) {
      DeltaCandidate* dc = &_data.listOfDeltaCandidates[idelta];
//-----------------------------------------------------------------------------
// don't extend candidates made out of one segment - but there is no such
// start from the first station to define limits
//-----------------------------------------------------------------------------
      int s1 = dc->fFirstStation;
      int s2 = dc->fLastStation-1;
      int last(-1);
//-----------------------------------------------------------------------------
// first check inside "holes"
//-----------------------------------------------------------------------------
      for (int i=s1; i<=s2; i++) {
        if (dc->seed[i] != nullptr) {
          last  = i;
          continue;
        }
//-----------------------------------------------------------------------------
// define expected T0 limits
//-----------------------------------------------------------------------------
        recoverStation(dc,last,i,1,0);
      }

      last  = dc->fFirstStation;
      for (int i=last-1; i>=0; i--) {
//-----------------------------------------------------------------------------
// skip empty stations
//-----------------------------------------------------------------------------
        if (dc->fFirstStation -i > _maxGap) break;
        recoverStation(dc,dc->fFirstStation,i,1,0);
      }

      last = dc->fLastStation;
      for (int i=last+1; i<kNStations; i++) {
//-----------------------------------------------------------------------------
// skip empty stations
//-----------------------------------------------------------------------------
        if (i-dc->fLastStation > _maxGap) break;
        recoverStation(dc,dc->fLastStation,i,1,0);
      }
    }

    return 0;
  }

//-----------------------------------------------------------------------------
// *FIXME* : need a formula for chi2, this simplification may not work well
// in the corners
//-----------------------------------------------------------------------------
  double DeltaFinder::seedDeltaChi2(DeltaSeed* Seed, DeltaCandidate* Delta) {

    int    nh   = Delta->NHits()+Seed->NHits();

    double nxym = (Delta->fSnxny+Seed->fSnxny)/nh;
    double nx2m = (Delta->fSnx2 +Seed->fSnx2 )/nh;
    double ny2m = (Delta->fSny2 +Seed->fSny2 )/nh;
    double nxrm = (Delta->fSnxnr+Seed->fSnxnr)/nh;
    double nyrm = (Delta->fSnynr+Seed->fSnynr)/nh;

    double d    = nx2m*ny2m-nxym*nxym;
    double xc   = (nyrm*nx2m-nxrm*nxym)/d;
    double yc   = (nyrm*nxym-nxrm*ny2m)/d;

    double dxs  = xc-Seed->Xc();
    double dys  = yc-Seed->Yc();

    double dxd  = xc-Delta->Xc();
    double dyd  = yc-Delta->Yc();

    double chi2 = (dxs*dxs+dys*dys+dxd*dxd+dyd*dyd)/(_maxDxy*_maxDxy);

    return chi2;
  }


//-----------------------------------------------------------------------------
// return 1 if a seed has been found , 0 otherwise
//-----------------------------------------------------------------------------
  int DeltaFinder::recoverSeed(DeltaCandidate* Delta, int LastStation, int Station) {
    // int rc(0);
                                        // predicted time range for this station
    float tdelta = Delta->T0(Station);
    // float xdelta = Delta->CofM.x();
    // float ydelta = Delta->CofM.y();

    float      chi2min(_maxChi2SeedDelta);
    DeltaSeed* closest_seed(nullptr);

    float dt   = _maxSeedDt + _maxDtDs*fabs(Station-LastStation);

    int nseeds = _data.NSeeds(Station);
    for (int i=0; i<nseeds; i++) {
      DeltaSeed* seed =  _data.deltaSeed(Station,i);
      if (seed->Good() == 0)                                          continue;
      if (seed->Used()     )                                          continue;
//-----------------------------------------------------------------------------
// one might need some safety here, but not the _maxDriftTime
//-----------------------------------------------------------------------------
      if (fabs(tdelta-seed->TMean()) > dt)                            continue;

      float chi2 = seedDeltaChi2(seed,Delta);

      if (chi2 < chi2min) {
                                        // new best seed
        closest_seed = seed;
        chi2min      = chi2;
      }
    }

    if (closest_seed) {
//-----------------------------------------------------------------------------
// the closest seed found, add it to the delta candidate and exit
// it is marked as associated with the delta candidate in DeltaCandidate::AddSeed
//-----------------------------------------------------------------------------
      Delta->AddSeed(closest_seed,Station);
      closest_seed->fChi2Delta = chi2min;
    }

    return (closest_seed != nullptr);
  }

//------------------------------------------------------------------------------
// try to recover hits of a 'Delta' candidate in a given 'Station'
// 'Delta' doesn't have hits in this station, check all hits here
// when predicting time, use the same value of Z for both layers of a given face
// return 1 if something has been found
//-----------------------------------------------------------------------------
  int DeltaFinder::recoverStation(DeltaCandidate* Delta, int LastStation, int Station, int UseUsedHits, int RecoverSeeds) {

                                        // predicted time range for this station
    float tdelta = Delta->T0(Station);
    float xdelta = Delta->CofM.x();
    float ydelta = Delta->CofM.y();
//-----------------------------------------------------------------------------
// first, loop over the existing seeds - need for the forward step
//-----------------------------------------------------------------------------
    if (RecoverSeeds) {
      int closest_seed_found = recoverSeed(Delta,LastStation,Station);
      if (closest_seed_found == 1) return 1;
    }
//-----------------------------------------------------------------------------
// no seeds found, look for single hits in the 'Station'
//-----------------------------------------------------------------------------
    DeltaSeed*  new_seed (nullptr);

    float dt_hit = _maxHitDt+_maxDtDs*fabs(Station-LastStation);

    for (int face=0; face<kNFaces; face++) {
      for (int ip=0; ip<kNPanelsPerFace; ip++) {
        PanelZ_t* panelz = &_data.oTracker[Station][face][ip];
        double dphi      = Delta->phi-panelz->phi;

        if (dphi < -M_PI) dphi += 2*M_PI;
        if (dphi >  M_PI) dphi -= 2*M_PI;
        if (fabs(dphi) > M_PI/3)                                      continue;

        if (tdelta-dt_hit > panelz->tmax)                             continue;
        if (tdelta+dt_hit < panelz->tmin)                             continue;
//-----------------------------------------------------------------------------
// panel and Delta overlap in phi and time, loop over hits
//-----------------------------------------------------------------------------
        int nhits = panelz->fHitData->size();
        for (int h=0; h<nhits; ++h) {
          HitData_t* hd = &(*panelz->fHitData)[h];
//-----------------------------------------------------------------------------
// don't skip hits already included into seeds - a two-hit stereo seed
// could be random
//-----------------------------------------------------------------------------
          if ((UseUsedHits == 0) and hd->Used())                      continue;
          float corr_time     = hd->fCorrTime;
//-----------------------------------------------------------------------------
// predicted time is the particle time, the drift time should be larger
//-----------------------------------------------------------------------------
          if (corr_time > tdelta+dt_hit)                              break;
          if (corr_time < tdelta-dt_hit)                              continue;

          const ComboHit*  ch = hd->fHit;
          double dx  = ch->pos().x()-xdelta;
          double dy  = ch->pos().y()-ydelta;
          double dw  = dx*panelz->wx+dy*panelz->wy; // distance along the wire
          double dxx = dx-panelz->wx*dw;
          double dyy = dy-panelz->wy*dw;

          double chi2_par  = (dw*dw)/(hd->fSigW2+_seedRes*_seedRes);
          double chi2_perp = (dxx*dxx+dyy*dyy)/_sigmaR2;
          double chi2      = chi2_par + chi2_perp;

          if (chi2 >= _maxChi2Radial)                                 continue;
          if (hd->Used()) {
//-----------------------------------------------------------------------------
// hit is a part of a seed. if the seed has 2 or less hits, don't check the chi2
// - that could be a random overlap
// if the seed has 3 or more hits, check the chi2
//-----------------------------------------------------------------------------
            int nh = hd->fSeed->NHits();
            if ((nh >= 3) and (chi2 > hd->fChi2Min))                   continue;
          }
//-----------------------------------------------------------------------------
// new hit needs to be added, create a special 1-hit seed for that
// in most cases, expect this seed not to have the second hit, but it may
// such a seed has its own CofM undefined
// ** FIXME ..in principle, at this point may want to check if the hit was used...
//-----------------------------------------------------------------------------
          if (new_seed == nullptr) {
            hd->fChi2Min = chi2;
            new_seed = _data.NewDeltaSeed(Station,face,hd,-1,nullptr);
          }
          else {
            if (face == new_seed->SFace(0)) {
//-----------------------------------------------------------------------------
// another close hit in the same panel, choose the best
//-----------------------------------------------------------------------------
              if (chi2 >= new_seed->HitData(face)->fChi2Min)          continue;
//-----------------------------------------------------------------------------
// new best hit in the same face
//-----------------------------------------------------------------------------
              hd->fChi2Min = chi2;
              new_seed->ReplaceFirstHit(hd);
            }
            else {
//-----------------------------------------------------------------------------
// more than one hit added in the hit pickup mode. The question is why a seed,
// constructed out of those two hits has not been added (or created)
//-----------------------------------------------------------------------------
              if (_printErrors) {
                printf("ERROR in DeltaFinder::recoverStation: ");
                printf("station=%2i - shouldn\'t be getting here, printout of new_seed and hd follows\n",Station);
                printf("chi2_par, chi2_perp, chi2: %8.2f %8.2f %8.2f\n",chi2_par, chi2_perp, chi2);

                printf("DELTA:\n");
                _data.printDeltaCandidate(Delta,"");
                printf("SEED:\n");
                _data.printDeltaSeed(new_seed,"");
                printf("HIT:\n");
                _data.printHitData  (hd      ,"");
              }

              new_seed->AddHit(hd,face);
            }

            if (corr_time < new_seed->fMinHitTime) new_seed->fMinHitTime = corr_time;
            if (corr_time > new_seed->fMaxHitTime) new_seed->fMaxHitTime = corr_time;
          }
        }
      }
      if (new_seed) new_seed->fFaceProcessed[face] = 1;
    }
//-----------------------------------------------------------------------------
// station is processed, see if anything has been found
// some parameters of seeds found in a recovery mode are not defined because
// there was no pre-seeding, for example
//-----------------------------------------------------------------------------
    int rc(0);
    if (new_seed) {
      int face0                       = new_seed->SFace(0);
      new_seed->HitData(face0)->fSeed = new_seed;

      Delta->AddSeed(new_seed,Station);
      rc = 1;
    }
                                        // return 1 if hits were added
    return rc;
  }

//-----------------------------------------------------------------------------
  void  DeltaFinder::runDeltaFinder() {

    for (int is=0; is<kNStations; is++) {
      for (int face=0; face<kNFaces; face++) {
        for (int ip=0; ip<kNPanelsPerFace; ip++) {
          PanelZ_t* pz = &_data.oTracker[is][face][ip];
          pz->fHitData->clear() ;
          pz->tmin =  1.e6;
          pz->tmax = -1.e6;
        }
      }
    }
    orderHits();
//-----------------------------------------------------------------------------
// loop over all stations and find delta seeds - 2-3-4 combo hit stubs
// a seed is always a stereo object
//-----------------------------------------------------------------------------
    findSeeds();
//-----------------------------------------------------------------------------
// connect seeds and create delta candidates
// at this stage, extend seeds to pick up single its in neighbor stations
// for single hits do not allo gaps
//-----------------------------------------------------------------------------
    connectSeeds();
//-----------------------------------------------------------------------------
// for existing delta candidates, pick up single gap-separated single hits
// no new candidates is created at this step
//-----------------------------------------------------------------------------
    recoverMissingHits();
//-----------------------------------------------------------------------------
// after recovery of missing hits, it is possible that some of delta candidates
// may need to be merged - try again
//-----------------------------------------------------------------------------
    mergeDeltaCandidates();
  }
}

//-----------------------------------------------------------------------------
// magic that makes this class a module.
//-----------------------------------------------------------------------------
DEFINE_ART_MODULE(mu2e::DeltaFinder)
//-----------------------------------------------------------------------------
// done
//-----------------------------------------------------------------------------
