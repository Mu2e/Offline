///////////////////////////////////////////////////////////////////////////////
// AgnosticHelixFinder
// M. Stortini, E. Martinez, N. Mazotov
///////////////////////////////////////////////////////////////////////////////

#include "Offline/CalPatRec/inc/AgnosticHelixFinder_types.hh"

#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Utilities/make_tool.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Sequence.h"

#include "Offline/BFieldGeom/inc/BFieldManager.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"

#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/HelixHit.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/RecoDataProducts/inc/StrawHitIndex.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/RecoDataProducts/inc/TrkFitDirection.hh"
#include "Offline/RecoDataProducts/inc/TrkFitFlag.hh"
#include "Offline/RecoDataProducts/inc/HelixRecoDir.hh"

#include "Offline/Mu2eUtilities/inc/LsqSums2.hh"
#include "Offline/Mu2eUtilities/inc/LsqSums4.hh"
#include "Offline/Mu2eUtilities/inc/polyAtan2.hh"
#include "Offline/Mu2eUtilities/inc/HelixTool.hh"
#include "Offline/Mu2eUtilities/inc/StopWatch.hh"

#include "Offline/Mu2eUtilities/inc/ModuleHistToolBase.hh"

#include "CLHEP/Units/PhysicalConstants.h"

namespace mu2e {

  using namespace AgnosticHelixFinderTypes;

  class AgnosticHelixFinder : public art::EDProducer {

  public:
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<int>             diagLevel              {Name("diagLevel"            ), Comment("turn tool on or off"         )  };
      fhicl::Atom<int>             doTiming               {Name("doTiming"             ), Comment("Analyze processing time"     ), 0};
      fhicl::Atom<art::InputTag>   chCollLabel            {Name("chCollLabel"          ), Comment("combo hit collection label"  )  };
      fhicl::Atom<art::InputTag>   tcCollLabel            {Name("tcCollLabel"          ), Comment("time cluster coll label"     )  };
      fhicl::Atom<art::InputTag>   ccCollLabel            {Name("ccCollLabel"          ), Comment("Calo Cluster coll label"     )  };
      fhicl::Atom<bool>            findMultipleHelices    {Name("findMultipleHelices"  ), Comment("allow more than one helix"   )  };
      fhicl::Atom<bool>            useStoppingTarget      {Name("useStoppingTarget"    ), Comment("allow in triplet candidates" )  };
      fhicl::Atom<int>             intenseEventThresh     {Name("intenseEventThresh"   ), Comment("# of clusters threshold"     )  };
      fhicl::Atom<int>             intenseClusterThresh   {Name("intenseClusterThresh" ), Comment("# of combo hits threshold"   )  };
      fhicl::Atom<bool>            doIsolationFlag        {Name("doIsolationFlag"      ), Comment("to filter out isolated hits" )  };
      fhicl::Atom<float>           isoRad                 {Name("isoRad"               ), Comment("for isolation cut"           )  };
      fhicl::Atom<int>             isoMinHitsNear         {Name("isoMinHitsNear"       ), Comment("#hits threshold for iso cut" )  };
      fhicl::Atom<bool>            doAverageFlag          {Name("doAverageFlag"        ), Comment("to average out hits or not"  )  };
      fhicl::Atom<bool>            doEDepFlag             {Name("doEDepFlag"           ), Comment("avoid protons in tripleting" )  };
      fhicl::Atom<float>           eDepFlagThresh         {Name("eDepFlagThresh"       ), Comment("threshold to set eDep flag"  )  };
      fhicl::Atom<float>           minDistCut             {Name("minDistCut"           ), Comment("for averaging out points"    )  };
      fhicl::Atom<float>           minTripletSeedZ        {Name("minTripletSeedZ"      ), Comment("minimum z for triplet seed"  )  };
      fhicl::Atom<float>           minTripletDz           {Name("minTripletDz"         ), Comment("min Z dist btwn 2 trip pts"  )  };
      fhicl::Atom<float>           maxTripletDz           {Name("maxTripletDz"         ), Comment("max Z dist btwn 2 trip pts"  )  };
      fhicl::Atom<float>           minTripletDist         {Name("minTripletDist"       ), Comment("min XY dist btwn 2 trip pts" )  };
      fhicl::Atom<float>           minTripletArea         {Name("minTripletArea"       ), Comment("triangle area of triplet"    )  };
      fhicl::Atom<float>           maxSeedCircleResidual  {Name("maxSeedCircleResidual"), Comment("add hits to triplet circle"  )  };
      fhicl::Atom<int>             minSeedCircleHits      {Name("minSeedCircleHits"    ), Comment("min hits to continue search" )  };
      fhicl::Atom<int>             maxSeedCircleHits      {Name("maxSeedCircleHits"    ), Comment("max hits to include in the seed" ), -1  };
      fhicl::Atom<float>           maxDphiDz              {Name("maxDphiDz"            ), Comment("used finding phi-z segment"  )  };
      fhicl::Atom<float>           maxSeedLineGap         {Name("maxSeedLineGap"       ), Comment("used finding phi-z segment"  )  };
      fhicl::Atom<float>           maxZWindow             {Name("maxZWindow"           ), Comment("used finding phi-z segment"  )  };
      fhicl::Atom<int>             minLineSegmentHits     {Name("minLineSegmentHits"   ), Comment("used in findSeedPhiLines()"  )  };
      fhicl::Atom<float>           segMultiplier          {Name("segMultiplier"        ), Comment("used in findSeedPhiLines()"  )  };
      fhicl::Atom<float>           maxSegmentChi2         {Name("maxSegmentChi2"       ), Comment("used in findSeedPhiLines()"  )  };
      fhicl::Atom<float>           max2PiAmbigResidual    {Name("max2PiAmbigResidual"  ), Comment("when 2pi resolving segment"  )  };
      fhicl::Atom<float>           maxPhiZResidual        {Name("maxPhiZResidual"      ), Comment("when refining phi-z line"    )  };
      fhicl::Atom<int>             minFinalSeedHits       {Name("minFinalSeedHits"     ), Comment("halt search if below thresh" )  };
      fhicl::Atom<float>           maxNHitsRatio          {Name("maxNHitsRatio"        ), Comment("max ratio of seed hits"      )  };
      fhicl::Atom<float>           maxCircleRecoverSigma  {Name("maxCircleRecoverSigma"), Comment("when doing final recovery"   )  };
      fhicl::Atom<float>           maxLineRecoverSigma    {Name("maxLineRecoverSigma"  ), Comment("when doing final recovery"   )  };
      fhicl::Atom<float>           caloClusterSigma       {Name("caloClusterSigma"     ), Comment("error assigned to calo clust")  };
      fhicl::Atom<int>             minNHelixStrawHits     {Name("minNHelixStrawHits"   ), Comment("straw hit save threshold"    )  };
      fhicl::Atom<int>             minNHelixComboHits     {Name("minNHelixComboHits"   ), Comment("combo hit save threshold"    )  };
      fhicl::Atom<float>           minHelixPerpMomentum   {Name("minHelixPerpMomentum" ), Comment("min pt of helix"             )  };
      fhicl::Atom<float>           maxHelixPerpMomentum   {Name("maxHelixPerpMomentum" ), Comment("max pt of helix"             )  };
      fhicl::Atom<float>           minHelixMomentum       {Name("minHelixMomentum"     ), Comment("min momentum of helix"       )  };
      fhicl::Atom<float>           maxHelixMomentum       {Name("maxHelixMomentum"     ), Comment("max momentum of helix"       )  };
      fhicl::Atom<float>           chi2LineSaveThresh     {Name("chi2LineSaveThresh"   ), Comment("max chi2Dof for line"        )  };
      fhicl::Atom<float>           maxEDepAvg             {Name("maxEDepAvg"           ), Comment("max avg edep of combohits"   )  };
      fhicl::Atom<float>           tzSlopeSigThresh       {Name("tzSlopeSigThresh"     ), Comment("direction ambiguous if below")  };
      fhicl::Sequence<std::string> validHelixDirections   {Name("validHelixDirections" ), Comment("only save desired directions")  };

      fhicl::Table<AgnosticHelixFinderTypes::Config> diagPlugin  {Name("diagPlugin"), Comment("diag plugin"                   )  };
    };

  private:
    //-----------------------------------------------------------------------------
    // tool on or off
    //-----------------------------------------------------------------------------
    int _diagLevel;
    int _doTiming;
    std::unique_ptr<StopWatch> _watch;

    //-----------------------------------------------------------------------------
    // event object labels
    //-----------------------------------------------------------------------------
    art::InputTag _chLabel;
    art::InputTag _tcLabel;
    art::InputTag _ccLabel;

    //-----------------------------------------------------------------------------
    // collections
    //-----------------------------------------------------------------------------
    const ComboHitCollection*      _chColl;
    const TimeClusterCollection*   _tcColl;
    const CaloClusterCollection*   _ccColl;

    //-----------------------------------------------------------------------------
    // helix search parameters
    //-----------------------------------------------------------------------------
    bool     _findMultipleHelices;
    bool     _useStoppingTarget;
    int      _intenseEventThresh;
    int      _intenseClusterThresh;
    bool     _doIsolationFlag;
    float    _isoRad;
    float    _isoRadSq;
    int      _isoMinHitsNear;
    bool     _doAverageFlag;
    bool     _doEDepFlag;
    float    _eDepFlagThresh;
    float    _minDistCut;
    float    _minDistCutSq;
    float    _minTripletSeedZ;
    float    _minTripletDz;
    float    _maxTripletDz;
    float    _minTripletDist;
    float    _minTripletDistSq;
    float    _minTripletArea;
    float    _maxSeedCircleResidual;
    float    _maxSeedCircleResidualSq;
    int      _minSeedCircleHits;
    int      _maxSeedCircleHits;
    float    _maxDphiDz;
    float    _maxSeedLineGap;
    float    _maxZWindow;
    int      _minLineSegmentHits;
    float    _segMultiplier;
    float    _segMultiplierSq;
    float    _maxSegmentChi2;
    float    _max2PiAmbigResidual;
    float    _maxPhiZResidual;
    int      _minFinalSeedHits;
    float    _maxNHitsRatio;
    float    _maxCircleRecoverSigma;
    float    _maxLineRecoverSigma;
    float    _caloClusterSigma;
    float    _caloClusterSigmaSq;
    int      _minNHelixStrawHits;
    int      _minNHelixComboHits;
    float    _minHelixPerpMomentum;
    float    _maxHelixPerpMomentum;
    float    _minHelixMomentum;
    float    _minHelixMomentumSq;
    float    _maxHelixMomentum;
    float    _maxHelixMomentumSq;
    float    _chi2LineSaveThresh;
    float    _maxEDepAvg;
    float    _tzSlopeSigThresh;
    std::vector<TrkFitDirection::FitDirection> _validHelixDirections;

    //-----------------------------------------------------------------------------
    // diagnostics
    //-----------------------------------------------------------------------------
    art::Handle<CaloClusterCollection>   _ccHandle;
    const mu2e::Calorimeter*             _calorimeter;
    const mu2e::Tracker*                 _tracker;
    art::Handle<TimeClusterCollection>   _tcCollH;

    //-----------------------------------------------------------------------------
    // stuff for doing helix search
    //-----------------------------------------------------------------------------
    std::vector<cHit>             _tcHits;
    XYZVectorF                    _stopTargPos;
    XYZVectorF                    _caloPos;
    ::LsqSums4                    _circleFitter;
    ::LsqSums2                    _lineFitter;
    std::vector<lineInfo>         _seedPhiLines;
    float                         _bz0;
    float                         _bz0Conv; // with conversion factor
    float                         _bz0ConvSq; // with conversion factor, squared
    float                         _caloZOffset;
    bool                          _intenseEvent;
    bool                          _intenseCluster;

    //-----------------------------------------------------------------------------
    // stuff for tool
    //-----------------------------------------------------------------------------
    diagInfo                              _diagInfo;
    std::unique_ptr<ModuleHistToolBase>   _hmanager;

    //-----------------------------------------------------------------------------
    // constants
    //-----------------------------------------------------------------------------
    static constexpr float _mmTconversion   = CLHEP::c_light/1000.0;
    static constexpr double        _twopi   = 2.*M_PI;

    //-----------------------------------------------------------------------------
    // functions
    //-----------------------------------------------------------------------------

  public:
    explicit AgnosticHelixFinder(const art::EDProducer::Table<Config>& config);
    virtual ~AgnosticHelixFinder();

    virtual void beginJob  ();
    virtual void beginRun  (art::Run&);
    virtual void produce   (art::Event& e);
    virtual void endJob    ();

    //-----------------------------------------------------------------------------
    // helper functions
    //-----------------------------------------------------------------------------
    bool         findData                  (const art::Event& e);
    const XYZVectorF& getPos               (size_t tcHitsIndex);
    void         computeCircleError2       (size_t tcHitsIndex, float xC, float yC);
    float        computeCircleResidual2    (size_t tcHitsIndex, float xC, float yC, float rC);
    void         computeHelixPhi           (size_t tcHitsIndex, float xC, float yC);
    float        computeHelixMomentum2     (float radius, float dphidz);
    float        computeHelixPerpMomentum  (float radius);
    void         tcHitsFill                (size_t tc);
    void         setFlags                  ();
    void         resetFlags                ();
    bool         findHelix                 (size_t tc, HelixSeedCollection& HSColl);
    bool         passesFlags               (size_t tcHitsIndex);
    void         setTripletI               (size_t tcHitsIndex, triplet& trip, LoopCondition& outcome);
    void         setTripletJ               (size_t tcHitsIndex, triplet& trip, LoopCondition& outcome);
    void         setTripletK               (size_t tcHitsIndex, triplet& trip, LoopCondition& outcome);
    void         initTriplet               (const triplet& trip, LoopCondition& outcome);
    void         initSeedCircle            (LoopCondition& outcome);
    void         initHelixPhi              ();
    void         findSeedPhiLines          (LoopCondition& outcome);
    float        underEstimateSlope        (float phi1, float phi1Err2, float phi2, float phi2Err2, float dz);
    void         resolve2PiAmbiguities     ();
    bool         refinePhiLine             (size_t lineIndex);
    void         initFinalSeed             (LoopCondition& outcome);
    bool         recoverPoints             ();
    void         checkHelixViability       (LoopCondition& outcome);
    void         saveHelix                 (size_t tc, HelixSeedCollection& HSColl);
    bool         validHelixDirection       (TrkFitDirection::FitDirection direction);

  };

  //-----------------------------------------------------------------------------
  // module constructor
  //-----------------------------------------------------------------------------
  AgnosticHelixFinder::AgnosticHelixFinder(const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config},
    _diagLevel                     (config().diagLevel()                             ),
    _doTiming                      (config().doTiming()                              ),
    _chLabel                       (config().chCollLabel()                           ),
    _tcLabel                       (config().tcCollLabel()                           ),
    _ccLabel                       (config().ccCollLabel()                           ),
    _findMultipleHelices           (config().findMultipleHelices()                   ),
    _useStoppingTarget             (config().useStoppingTarget()                     ),
    _intenseEventThresh            (config().intenseEventThresh()                    ),
    _intenseClusterThresh          (config().intenseClusterThresh()                  ),
    _doIsolationFlag               (config().doIsolationFlag()                       ),
    _isoRad                        (config().isoRad()                                ),
    _isoMinHitsNear                (config().isoMinHitsNear()                        ),
    _doAverageFlag                 (config().doAverageFlag()                         ),
    _doEDepFlag                    (config().doEDepFlag()                            ),
    _eDepFlagThresh                (config().eDepFlagThresh()                        ),
    _minDistCut                    (config().minDistCut()                            ),
    _minTripletSeedZ               (config().minTripletSeedZ()                       ),
    _minTripletDz                  (config().minTripletDz()                          ),
    _maxTripletDz                  (config().maxTripletDz()                          ),
    _minTripletDist                (config().minTripletDist()                        ),
    _minTripletArea                (config().minTripletArea()                        ),
    _maxSeedCircleResidual         (config().maxSeedCircleResidual()                 ),
    _minSeedCircleHits             (config().minSeedCircleHits()                     ),
    _maxSeedCircleHits             (config().maxSeedCircleHits()                     ),
    _maxDphiDz                     (config().maxDphiDz()                             ),
    _maxSeedLineGap                (config().maxSeedLineGap()                        ),
    _maxZWindow                    (config().maxZWindow()                            ),
    _minLineSegmentHits            (config().minLineSegmentHits()                    ),
    _segMultiplier                 (config().segMultiplier()                         ),
    _maxSegmentChi2                (config().maxSegmentChi2()                        ),
    _max2PiAmbigResidual           (config().max2PiAmbigResidual()                   ),
    _maxPhiZResidual               (config().maxPhiZResidual()                       ),
    _minFinalSeedHits              (config().minFinalSeedHits()                      ),
    _maxNHitsRatio                 (config().maxNHitsRatio()                         ),
    _maxCircleRecoverSigma         (config().maxCircleRecoverSigma()                 ),
    _maxLineRecoverSigma           (config().maxLineRecoverSigma()                   ),
    _caloClusterSigma              (config().caloClusterSigma()                      ),
    _minNHelixStrawHits            (config().minNHelixStrawHits()                    ),
    _minNHelixComboHits            (config().minNHelixComboHits()                    ),
    _minHelixPerpMomentum          (config().minHelixPerpMomentum()                  ),
    _maxHelixPerpMomentum          (config().maxHelixPerpMomentum()                  ),
    _minHelixMomentum              (config().minHelixMomentum()                      ),
    _maxHelixMomentum              (config().maxHelixMomentum()                      ),
    _chi2LineSaveThresh            (config().chi2LineSaveThresh()                    ),
    _maxEDepAvg                    (config().maxEDepAvg()                            ),
    _tzSlopeSigThresh              (config().tzSlopeSigThresh()                      )
  {
    // to cache the evaluations
    _minTripletDistSq = _minTripletDist * _minTripletDist;
    _maxSeedCircleResidualSq = _maxSeedCircleResidual * _maxSeedCircleResidual;
    _caloClusterSigmaSq = _caloClusterSigma * _caloClusterSigma;
    _isoRadSq = _isoRad * _isoRad;
    _minDistCutSq = _minDistCut * _minDistCut;
    _minHelixMomentumSq = _minHelixMomentum * _minHelixMomentum;
    _maxHelixMomentumSq = _maxHelixMomentum * _maxHelixMomentum;
    _segMultiplierSq = _segMultiplier * _segMultiplier;

    // reserve a reasonable amount of space for future hit collections
    _tcHits.reserve(200);

      // convert the helix direction names into enums
    for(auto helix_dir : config().validHelixDirections()) {
      _validHelixDirections.push_back(TrkFitDirection::fitDirectionFromName(helix_dir));
    }
    consumes<ComboHitCollection>     (_chLabel);
    consumes<TimeClusterCollection>  (_tcLabel);
    consumes<CaloClusterCollection>  (_ccLabel);
    produces<HelixSeedCollection>    ();

    if (_diagLevel > 0) _hmanager = art::make_tool<ModuleHistToolBase>(config().diagPlugin, "diagPlugin");
    else _hmanager = std::make_unique<ModuleHistToolBase>();

    // For timing analysis
    if(_doTiming) _watch = std::make_unique<StopWatch>();

    if (_useStoppingTarget == true) { _stopTargPos.SetCoordinates(0.0, 0.0, std::numeric_limits<float>::max()); }

  }

  //-----------------------------------------------------------------------------
  // destructor
  //-----------------------------------------------------------------------------
  AgnosticHelixFinder::~AgnosticHelixFinder() {}

  //-----------------------------------------------------------------------------
  // beginJob
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::beginJob() {
    if (_diagLevel > 0) {
      art::ServiceHandle<art::TFileService> tfs;
      _hmanager->bookHistograms(tfs);
    }
    // calibrate the timing speed
    if(_doTiming) {
      printf("[AgnosticHelixFinder::%s] Calibrating timing...\n", __func__);
      _watch->Calibrate();
      printf("[AgnosticHelixFinder::%s] Calibrating time = %.3g us\n", __func__, _watch->Calibration());
    }
  }

  //-----------------------------------------------------------------------------
  // beginRun
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::beginRun(art::Run&) {

    GeomHandle<mu2e::Calorimeter> ch;
    _calorimeter = ch.get();

    GeomHandle<mu2e::Tracker> th;
    _tracker = th.get();

    GeomHandle<BFieldManager> bfmgr;
    GeomHandle<DetectorSystem> det;
    CLHEP::Hep3Vector vpoint_mu2e = det->toMu2e(CLHEP::Hep3Vector(0.0, 0.0, 0.0));
    _bz0 = bfmgr->getBField(vpoint_mu2e).z();
    _bz0Conv = _bz0 * _mmTconversion;
    _bz0ConvSq = _bz0Conv * _bz0Conv;

    // Offset for calo cluster z positions
    double offset = _calorimeter->caloInfo().getDouble("diskCaseZLength");
    offset += _calorimeter->caloInfo().getDouble("BPPipeZOffset");
    offset += _calorimeter->caloInfo().getDouble("BPHoleZLength");
    offset += _calorimeter->caloInfo().getDouble("FEEZLength");
    offset /= 2.0;
    _caloZOffset = offset;
  }

  //-----------------------------------------------------------------------------
  // find input data
  //-----------------------------------------------------------------------------
  bool AgnosticHelixFinder::findData(const art::Event& evt) {
    if(_doTiming) _watch->Increment(__func__);

    auto chCollH = evt.getValidHandle<ComboHitCollection>(_chLabel);
    _chColl = chCollH.product();

    auto tcCollH = evt.getValidHandle<TimeClusterCollection>(_tcLabel);
    _tcColl = tcCollH.product();
     // save for art::Ptr creation
    evt.getByLabel(_tcLabel, _tcCollH);

    // Calorimeter cluster collection is no longer used
    // if (evt.getByLabel(_ccLabel, _ccHandle)) {
    //   _ccColl = _ccHandle.product();
    // } else { _ccColl = nullptr; }


    // prepare diagnostic tool data members
    if (_diagLevel > 0) {
      _diagInfo.nHelices = 0;
      _diagInfo.timeClusterData.clear();
      _diagInfo.helixSeedData.clear();
      _diagInfo.lineSegmentData.clear();
      _diagInfo.hseed = nullptr;
      _diagInfo.event = &evt;
      _diagInfo.tcHits = &_tcHits;
      _diagInfo.seedPhiLines = &_seedPhiLines;
      _diagInfo.circleFitter = &_circleFitter;
      _diagInfo.lineFitter = &_lineFitter;
      _diagInfo.chColl = _chColl;
      _diagInfo.tcColl = _tcColl;
      _diagInfo.diagLevel = _diagLevel;
      _diagInfo.bz0 = _bz0;
      _hmanager->fillHistograms(&_diagInfo, DIAG::kBegin);
    }
    if(_doTiming) _watch->StopTime(__func__);

    return (_tcColl) && (_chColl); // require a time cluster collection and a combo hit collection
  }

  //-----------------------------------------------------------------------------
  // event entry point
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::produce(art::Event& event) {
    if(_doTiming) _watch->Increment(__func__);

    // Prepare the helix collection
    std::unique_ptr<HelixSeedCollection> hsColl(new HelixSeedCollection);

    // get the necessary data to do helix search and then search for helices
    if (findData(event)) {

      // flag whether event is intense or not
      _intenseEvent = (int) _tcColl->size() > _intenseEventThresh;
      if(_diagLevel > 1)
        printf("[AgnosticHelixFinder::%s] Event %i:%i:%i: N(time clusters) = %zu\n", __func__, event.run(), event.subRun(), event.event(), _tcColl->size());

      for (size_t i = 0; i < _tcColl->size(); i++) {
        if(_doTiming) _watch->Increment("per-time cluster");
        // check to see if cluster is a busy one
        // if cluster is busy or event is intense then only process the time cluster if it has calo cluster
        const auto& tc = _tcColl->at(i);
        _intenseCluster = (int) tc.nhits() > _intenseClusterThresh;
        if (_intenseEvent || _intenseCluster) {
          if (!tc.hasCaloCluster()) {
            if(_diagLevel > 1)
              printf("  --> Skipping time cluster (%zu hits) without a calo cluster\n", tc.nhits());
            continue;
          }
        }
        const int nHelicesInitial = _diagInfo.nHelices; // Only valid if diagLevel > 0, for diagnostic tracking
        tcHitsFill(i); // Initialize the list of hits in the time cluster
        while(findHelix(i, *hsColl) && _findMultipleHelices); // Exit the search if no helix is found or after finding a helix if not configured for multi-helix reco

        if (_diagLevel > 0) {
          tcInfo timeClusterInfo;
          timeClusterInfo.nHelices = _diagInfo.nHelices - nHelicesInitial;
          timeClusterInfo.nComboHits = _tcColl->at(i).nhits();
          timeClusterInfo.nStrawHits = _tcColl->at(i).nStrawHits();
          _diagInfo.timeClusterData.push_back(timeClusterInfo);
          _hmanager->fillHistograms(&_diagInfo, DIAG::kTimeCluster);
          _diagInfo.helixSeedData.clear(); // clear for the next time cluster
        }
      }
    } else if(_diagLevel > 0) {
      printf("[AgnosticHelixFinder::%s] %i:%i:%i Event data not found!\n", __func__, event.run(), event.subRun(), event.event());
    }
    if(_doTiming) _watch->StopTime("per-time cluster");

    // fill necessary data members for diagnostic tool
    if (_diagLevel > 0) {
      _diagInfo.nTimeClusters = _tcColl->size();
      for(auto& hseed : *hsColl) _diagInfo.helixSeedData.push_back(hsInfo(_bz0, hseed));
      _hmanager->fillHistograms(&_diagInfo, DIAG::kEnd);
    }

    // put helix seed collection into the event record
    if(_doTiming) _watch->Increment("output");
    event.put(std::move(hsColl));
    if(_doTiming) _watch->StopTime ("output");

    if(_doTiming) _watch->StopTime(__func__);
  }

  //-----------------------------------------------------------------------------
  // endJob
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::endJob() {
    if(_doTiming) {
      for(int itest = 0; itest < 100000; ++itest) _watch->Increment("AAA-TimeTest"); // put at the top a reference for timing impacts
      _watch->StopTime("AAA-TimeTest");
      std::cout << "[AgnosticHelixFinder::" << __func__ << "::" <<
        moduleDescription().moduleLabel() << "] Timing:\n" << *_watch;
    }
  }

  //-----------------------------------------------------------------------------
  // get pointer to position of hit given some index in _tcHits
  //-----------------------------------------------------------------------------
  inline const XYZVectorF& AgnosticHelixFinder::getPos(size_t tcHitsIndex) {
    return _tcHits[tcHitsIndex].pos;
    // const int hitIndice = _tcHits[tcHitsIndex].hitIndice;
    // if (hitIndice == HitType::STOPPINGTARGET) { return _stopTargPos; }
    // if (hitIndice == HitType::CALOCLUSTER) { return _caloPos; }
    // return _chColl->at(hitIndice).pos();
  }

  //-----------------------------------------------------------------------------
  // updates circleError of hit given some circle parameters
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::computeCircleError2(size_t tcHitsIndex, float xC, float yC) {

    auto& tcHit = _tcHits[tcHitsIndex];
    int hitIndice = tcHit.hitIndice;

    if (hitIndice == HitType::CALOCLUSTER) {
      tcHit.circleError2 = _caloClusterSigmaSq;
    } else if(hitIndice >= 0) {
      const auto& hit = *(tcHit.hit);
      float transVar = hit.transVar();
      const auto& pos = hit.pos();
      float x = pos.x();
      float y = pos.y();
      float dx = x - xC;
      float dy = y - yC;
      const auto& vdir = hit.vDir();
      float dxn = dx * vdir.x() + dy *vdir.y();
      float costh2 = dxn * dxn / (dx * dx + dy * dy);
      float sinth2 = 1.f - costh2;
      tcHit.circleError2 = hit.wireVar()*sinth2 + transVar*costh2;
    }
  }

  //-----------------------------------------------------------------------------
  // compute residual between point an circle
  //-----------------------------------------------------------------------------
  float AgnosticHelixFinder::computeCircleResidual2(size_t tcHitsIndex, float xC, float yC, float rC) {

    const auto& pos = getPos(tcHitsIndex);
    const float x = pos.x() - xC;
    const float y = pos.y() - yC;
    const float deltaDistance = rC - std::sqrt(x*x + y*y); // sign lost since it's squared
    const float circleSigma2 = _tcHits[tcHitsIndex].circleError2;

    return deltaDistance * deltaDistance / circleSigma2;
  }

  //-----------------------------------------------------------------------------
  // compute phi relative to helix center, and set helixPhiError
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::computeHelixPhi(size_t tcHitsIndex, float xC, float yC) {

    auto& tcHit = _tcHits[tcHitsIndex];
    const int hitIndice = tcHit.hitIndice;
    const auto& pos = tcHit.pos;

    const float X = pos.x() - xC;
    const float Y = pos.y() - yC;
    const float Rsq = X * X + Y * Y;
    tcHit.helixPhi = polyAtan2(Y, X);
    if (tcHit.helixPhi < 0.) {
      tcHit.helixPhi += _twopi;
    }
    // find phi error and initialize it
    // find phi error by projecting errors onto vector tangent to circle
    float deltaS2(0.f);
    if (hitIndice == HitType::CALOCLUSTER) {
      deltaS2 = _caloClusterSigmaSq;
    } else {
      const auto& hit = *(tcHit.hit);
      const float R = std::sqrt(Rsq);
      const auto& udir = hit.uDir();
      const float ux = udir.x();
      const float uy = udir.y();
      float tanVecX = Y / R;
      float tanVecY = -X / R;
      float wireErr = hit.wireRes();
      float wireVecX = ux;
      float wireVecY = uy;
      float projWireErr = wireErr * (wireVecX * tanVecX + wireVecY * tanVecY);
      float transErr = hit.transRes();
      float transVecX = uy;
      float transVecY = -ux;
      float projTransErr = transErr * (transVecX * tanVecX + transVecY * tanVecY);
      deltaS2 = projWireErr * projWireErr + projTransErr * projTransErr;
    }
    tcHit.helixPhiError2 = deltaS2 / Rsq;
  }

  //-----------------------------------------------------------------------------
  // function to compute helix momentum^2 given some circle radius and line slope
  //-----------------------------------------------------------------------------
  inline float AgnosticHelixFinder::computeHelixMomentum2(float radius, float dphidz) {

    float lambda = 1.f / dphidz;

    return (radius * radius + lambda * lambda) * _bz0ConvSq;
  }

  //-----------------------------------------------------------------------------
  // function to compute helix transverse momentum given circle radius and b-field
  //-----------------------------------------------------------------------------
  inline float AgnosticHelixFinder::computeHelixPerpMomentum(float radius) {

    return radius * _bz0Conv;
  }

  //-----------------------------------------------------------------------------
  // fill vector with hits to search for helix
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::tcHitsFill(size_t itc) {
    if(_doTiming) _watch->SetTime(__func__);
    _tcHits.clear(); // clear the initial list

    const auto& tc = _tcColl->at(itc);
    size_t sortStartIndex = 0;

    // push back stopping target if it is to be used
    if (_useStoppingTarget == true) {
      cHit hit;
      hit.hitIndice = HitType::STOPPINGTARGET;
      hit.pos = _stopTargPos;
      _tcHits.push_back(hit);
      sortStartIndex++;
    }

    // push back calo cluster if it exists in time cluster
    if (tc.hasCaloCluster()) {
      const art::Ptr<CaloCluster> cl = tc.caloCluster();
      if (cl.isNonnull()) {
        cHit hit;
        hit.hitIndice = HitType::CALOCLUSTER;
        CLHEP::Hep3Vector gpos = _calorimeter->geomUtil().diskToMu2e(cl->diskID(), cl->cog3Vector());
        CLHEP::Hep3Vector tpos = _calorimeter->geomUtil().mu2eToTracker(gpos);
        _caloPos.SetCoordinates(tpos.x(), tpos.y(), tpos.z() - _caloZOffset);
        hit.pos = _caloPos;
        _tcHits.push_back(hit);
        sortStartIndex++;
      }
    }

    // fill hits from time cluster
    for (size_t i = 0; i < tc._strawHitIdxs.size(); i++) {
      cHit hit;
      hit.hitIndice = tc._strawHitIdxs[i];
      hit.hit = &(_chColl->at(hit.hitIndice));
      hit.pos = hit.hit->pos();
      _tcHits.push_back(hit);
    }

    // order from largest z to smallest z (leave the stopping target and calo cluster hits at the front)
    std::sort(_tcHits.begin() + sortStartIndex, _tcHits.end(), [&](const cHit& a, const cHit& b) {
      return a.pos.z() > b.pos.z();
    });

    // For diagnostics
    if(_diagLevel > 0) {
      _diagInfo.targPos = _stopTargPos;
      _diagInfo.caloPos = _caloPos;
      _diagInfo.tc = &tc;
    }
    if(_doTiming) _watch->Increment(__func__);
  }

  //-----------------------------------------------------------------------------
  // logic for setting certain flags on hits
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::setFlags() {
    if(!_doIsolationFlag && !_doAverageFlag && !_doEDepFlag) return; // nothing to flag

    // do isolation, average, and eDepFlag flagging
    const size_t ntcHits = _tcHits.size();
    for (size_t i = 0; i < ntcHits; i++) {
      auto& tcHit_i = _tcHits[i];
      const int index_i = tcHit_i.hitIndice;
      if (tcHit_i.inHelix
          || index_i == HitType::STOPPINGTARGET
          || index_i == HitType::CALOCLUSTER) { continue; }

      // test for high energy deposition
      tcHit_i.highEDep = _doEDepFlag && tcHit_i.hit->energyDep() > _eDepFlagThresh;

      // check for isolated hits or average out nearby hits
      if (_doIsolationFlag == true || _doAverageFlag == true) {
        int nHitsNear = 0;
        const auto& seedPos = getPos(i);
        tcHit_i.isolated = _doIsolationFlag;
        for (size_t j = 0; j < ntcHits; j++) {
          if(j == i) continue; // don't check against itself
          if(!tcHit_i.isolated) { // passed isolation check already
            if(!_doAverageFlag || tcHit_i.averagedOut) break; // nothing left to check
            if(j < i) {  // skip ahead since early js would be checked against this hit already for averaging out
              j = i;
              continue;
            }
          }
          auto& tcHit_j = _tcHits[j];
          if(!tcHit_i.isolated && tcHit_j.averagedOut) continue; // nothing to do

          const int index_j = tcHit_j.hitIndice;
          if (tcHit_j.inHelix
              || index_j == HitType::STOPPINGTARGET
              || index_j == HitType::CALOCLUSTER) { continue; }
          const auto& testPos = getPos(j);
          const float perp2 = (std::pow(seedPos.x() - testPos.x(),2) + // compute by hand to avoid object creation
                               std::pow(seedPos.y() - testPos.y(),2));
          // do isolation flagging unless already identified as non-isolated
          if (tcHit_i.isolated) {
            if (perp2 < _isoRadSq) { nHitsNear++; } // if near the hit i, increment the counter
            tcHit_i.isolated = nHitsNear < _isoMinHitsNear;
          }
          // do averaging out
          if (_doAverageFlag == true) {
            if (tcHit_i.averagedOut == true) { continue; }
            if (tcHit_j.averagedOut == true) { continue; }
            if (perp2 <= _minDistCutSq) { tcHit_j.averagedOut = true; }
          }
        }
      }
    }
  }

  //-----------------------------------------------------------------------------
  // logic for setting certain flags on hits
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::resetFlags() {

    // reset flags that need to be reset before next helix search
    for (size_t i = 0; i < _tcHits.size(); i++) {
      if (_tcHits[i].inHelix == true || _tcHits[i].hitIndice == HitType::STOPPINGTARGET) { continue; }
      _tcHits[i].averagedOut = false;
      _tcHits[i].notOnLine = true;
      _tcHits[i].notOnSegment = true;
    }
  }

  //-----------------------------------------------------------------------------
  // logic to find helix
  //-----------------------------------------------------------------------------
  bool AgnosticHelixFinder::findHelix(size_t tc, HelixSeedCollection& HSColl) {
    if(_doTiming) _watch->SetTime(__func__);

    // clear line and circle fitters
    _circleFitter.clear();
    _lineFitter.clear();

    // set flags before starting search so that we know what hits to use for tripletting
    if (_doIsolationFlag == true || _doAverageFlag == true) { setFlags(); }

    // now we loop over triplets
    LoopCondition loopCondition;
    triplet tripletInfo;
    for (size_t i = 0; i < _tcHits.size() - 2; i++) {
      if(_doTiming > 1) _watch->Increment("triplet-i");
      if(_doTiming > 2) _watch->StopTime ("triplet-j");
      bool uselessSeed = true;
      setTripletI(i, tripletInfo, loopCondition);
      if (loopCondition == CONTINUE) { continue; }
      if (loopCondition == BREAK) { break; }
      for (size_t j = i + 1; j < _tcHits.size() - 1; j++) {
        if(_doTiming > 2) _watch->Increment("triplet-j");
        if(_doTiming > 3) _watch->StopTime ("triplet-k");
        setTripletJ(j, tripletInfo, loopCondition);
        if (loopCondition == CONTINUE) { continue; }
        if (loopCondition == BREAK) { break; }
        for (size_t k = j + 1; k < _tcHits.size(); k++) {
          if(_doTiming > 3) _watch->Increment("triplet-k");
          setTripletK(k, tripletInfo, loopCondition);
          if (loopCondition == CONTINUE) { continue; }
          if (loopCondition == BREAK) { break; }
          // clear fitters for new triplet then initialize triplet
          _circleFitter.clear();
          _lineFitter.clear();
          initTriplet(tripletInfo, loopCondition);
          // now initialize seed circle if triplet circle passed condition check
          if (loopCondition == CONTINUE) { continue; }
          initSeedCircle(loopCondition);
          if (loopCondition == CONTINUE) { continue; }
          uselessSeed = false;
          initHelixPhi();
          findSeedPhiLines(loopCondition);
          if (loopCondition == CONTINUE) { continue; }
          resolve2PiAmbiguities();
          // refine the seed phi lines by removing the worst hits
          for (size_t ii = 0; ii < _seedPhiLines.size(); ii++) {
            if ((int)_seedPhiLines[ii].tcHitsIndices.size() < _minFinalSeedHits) { continue; }
            while (refinePhiLine(ii)); // while refinements are still made, keep going
          }
          initFinalSeed(loopCondition);
          if (loopCondition == CONTINUE) { continue; }
          while (recoverPoints()); // while points are recovered, keep trying to recover
          checkHelixViability(loopCondition);
          if (loopCondition == CONTINUE) { continue; }
          // before saving helix we make sure it has enough hits
          int nStrawHitsInHelix = 0;
          int nComboHitsInHelix = 0;
          int nStrawHitsInTimeCluster = 0;
          int nComboHitsInTimeCluster = 0;
          // compute number of usable hits in time cluster, and number of hits in candidate helix
          for (size_t q = 0; q < _tcHits.size(); q++) {
            if (_tcHits[q].inHelix == true || _tcHits[q].hitIndice == HitType::STOPPINGTARGET
                || _tcHits[q].hitIndice == HitType::CALOCLUSTER) { continue; }
            const auto& tcHit = _tcHits[q];
            const int nStrawHitsInHit = tcHit.hit->nStrawHits();
            nStrawHitsInTimeCluster += nStrawHitsInHit;
            ++nComboHitsInTimeCluster;
            if (tcHit.used == false) { continue; }
            nStrawHitsInHelix += nStrawHitsInHit;
            ++nComboHitsInHelix;
          }
          if (nStrawHitsInHelix >= _minNHelixStrawHits && nComboHitsInHelix >= _minNHelixComboHits) {
            saveHelix(tc, HSColl);
            if (_diagLevel > 0) { _diagInfo.nHelices++; }
            // we only want to search for another helix if we have enough remaining hits after saving helix
            const int remainingStrawHits = nStrawHitsInTimeCluster - nStrawHitsInHelix;
            const int remainingComboHits = nComboHitsInTimeCluster - nComboHitsInHelix;
            const bool findAnotherHelix = remainingStrawHits >= _minNHelixStrawHits && remainingComboHits >= _minNHelixComboHits;
            if(findAnotherHelix && _doAverageFlag) resetFlags(); // reset flags if needed

            // Helix found, so we can exit this search
            if(_doTiming > 0) _watch->StopTime(__func__);
            if(_doTiming > 1) _watch->StopTime("triplet-i");
            if(_doTiming > 2) _watch->StopTime("triplet-j");
            if(_doTiming > 3) _watch->StopTime("triplet-k");
            return findAnotherHelix;
          } else if(_diagLevel > 4) {
            printf("Helix hits: Failed with %i straw hits (%i combo hits)\n", nStrawHitsInHelix, nComboHitsInHelix);
          }
        } // end triplet k loop
      } // end triplet j loop
      _tcHits[i].uselessTripletSeed = uselessSeed;
    } // end triplet i loop

    if(_doTiming > 0) _watch->StopTime(__func__);
    if(_doTiming > 1) _watch->StopTime("triplet-i");
    if(_doTiming > 2) _watch->StopTime("triplet-j");
    if(_doTiming > 3) _watch->StopTime("triplet-k");

    // No helix found, no need to continue searching in this time cluster
    return false;
  }

  //-----------------------------------------------------------------------------
  // check flags to see if point is good for triplet-ing with
  //-----------------------------------------------------------------------------
  bool AgnosticHelixFinder::passesFlags(size_t tcHitsIndex) {

    if (_tcHits[tcHitsIndex].inHelix == true) { return false; }
    if (_tcHits[tcHitsIndex].uselessTripletSeed == true) { return false; }
    if (_doIsolationFlag == true && _tcHits[tcHitsIndex].isolated == true) { return false; }
    if (_doAverageFlag == true && _tcHits[tcHitsIndex].averagedOut == true) { return false; }
    if (_doEDepFlag == true && _tcHits[tcHitsIndex].highEDep == true) { return false; }


    return true;
  }

  //-----------------------------------------------------------------------------
  // set ith point of triplet, return outcome value which can be used to direct for loops
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::setTripletI(size_t tcHitsIndex, triplet& trip, LoopCondition& outcome) {

    // don't clone the data until necessary
    const auto& pos = getPos(tcHitsIndex);
    const int index = _tcHits[tcHitsIndex].hitIndice;

    // check if point should break for loop
    if (pos.z() < _minTripletSeedZ) {
      outcome = BREAK;
      return;
    }
    if (_intenseEvent == true || _intenseCluster == true) {
      if (index != HitType::CALOCLUSTER && index != HitType::STOPPINGTARGET) {
        outcome = BREAK;
        return;
      }
    }

    // check if point should be continued on
    if (!passesFlags(tcHitsIndex)) {
      outcome = CONTINUE;
      return;
    }

    // set info for the point
    trip.i.pos = pos;
    trip.i.hitIndice = index;

    outcome = GOOD;
  }

  //-----------------------------------------------------------------------------
  // set jth point of triplet, return outcome value which can be used to direct for loops
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::setTripletJ(size_t tcHitsIndex, triplet& trip, LoopCondition& outcome) {
    // don't clone the data until necessary
    const auto& pos = getPos(tcHitsIndex);
    const int index = _tcHits[tcHitsIndex].hitIndice;

    const float dz12 = trip.i.pos.z() - pos.z();

    // check if point should break for loop
    if (trip.i.hitIndice != HitType::CALOCLUSTER && trip.i.hitIndice != HitType::STOPPINGTARGET){
        if (dz12 > _maxTripletDz) {
          outcome = BREAK;
          return;
        }
    }
    if (_intenseEvent == true || _intenseCluster == true) {
      if (trip.i.hitIndice == HitType::STOPPINGTARGET && index != HitType::CALOCLUSTER) {
          outcome = BREAK;
          return;
      }
    }

    // check if point should be continued on
    if (!passesFlags(tcHitsIndex) || dz12 < _minTripletDz ||
        (std::pow(trip.i.pos.x()-pos.x(),2) + std::pow(trip.i.pos.y() - pos.y(), 2)) < _minTripletDistSq) {
        // (trip.i.pos - pos).Perp2() < _minTripletDistSq) {
      outcome = CONTINUE;
      return;
    }

    // set info for the point
    trip.j.pos = pos;
    trip.j.hitIndice = index;

    outcome = GOOD;
  }

  //-----------------------------------------------------------------------------
  // set kth point of triplet, return outcome value which can be used to direct for loops
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::setTripletK(size_t tcHitsIndex, triplet& trip, LoopCondition& outcome) {

    // don't clone the data until necessary
    const auto& pos = getPos(tcHitsIndex);
    const int index = _tcHits[tcHitsIndex].hitIndice;

    const float dz23 = trip.j.pos.z() - pos.z();

    // check if point should break for loop
    if (trip.j.hitIndice != HitType::CALOCLUSTER && trip.j.hitIndice != HitType::STOPPINGTARGET){
        if(dz23 > _maxTripletDz) {
          outcome = BREAK;
          return;
        }
    }

    // check if point should be continued on
    if (!passesFlags(tcHitsIndex) || dz23 < _minTripletDz ||
        (trip.i.pos-pos).Perp2() < _minTripletDistSq ||
        (trip.j.pos-pos).Perp2() < _minTripletDistSq) {
      outcome = CONTINUE;
      return;
    }

    // set the info for point
    trip.k.pos = pos;
    trip.k.hitIndice = index;

    outcome = GOOD;
  }

  //-----------------------------------------------------------------------------
  // finding circle from triplet
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::initTriplet(const triplet& trip, LoopCondition& outcome) {
    if(_doTiming > 4) _watch->SetTime(__func__);
    _circleFitter.addPoint(trip.i.pos.x(), trip.i.pos.y());
    _circleFitter.addPoint(trip.j.pos.x(), trip.j.pos.y());
    _circleFitter.addPoint(trip.k.pos.x(), trip.k.pos.y());

    // check if circle is valid for search or if we should continue to next triplet
    float radius = _circleFitter.radius();
    float pt = computeHelixPerpMomentum(radius);
    if (pt < _minHelixPerpMomentum || pt > _maxHelixPerpMomentum) { outcome = CONTINUE; }
    else { outcome = GOOD; }

    if(_diagLevel > 0) { // diagnose triplets
      _diagInfo.tripInfo.trip = trip;
      _diagInfo.tripInfo.radius = radius;
      _diagInfo.tripInfo.xC = _circleFitter.x0();
      _diagInfo.tripInfo.yC = _circleFitter.y0();
      _diagInfo.loopCondition = outcome;

      _hmanager->fillHistograms(&_diagInfo, DIAG::kTriplet);
    }
    if(_doTiming > 4) _watch->StopTime(__func__);

  }

  //-----------------------------------------------------------------------------
  // start with initial seed circle
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::initSeedCircle(LoopCondition& outcome) {
    if(_doTiming > 4) _watch->SetTime(__func__);

    // get triplet circle parameters then clear fitter
    float xC = _circleFitter.x0();
    float yC = _circleFitter.y0();
    float rC = _circleFitter.radius();
    _circleFitter.clear();

    // project error bars onto the triplet circle found and add to fitter those within defined max
    // residual
    for (size_t i = 0; i < _tcHits.size(); i++) {
      auto& hit = _tcHits[i];
      if (hit.inHelix) { continue; }
      if (hit.hitIndice == HitType::STOPPINGTARGET) {
        hit.used = false;
        continue;
      }
      if(_doTiming > 5) _watch->SetTime("initSeedCircle-error");
      computeCircleError2(i, xC, yC);
      if(_doTiming > 5) _watch->StopTime("initSeedCircle-error");
      if ((_maxSeedCircleHits < 0 || _circleFitter.qn() < _maxSeedCircleHits) &&
          computeCircleResidual2(i, xC, yC, rC) < _maxSeedCircleResidualSq) {
        const float wP = 1.f / (hit.circleError2);
        const auto& pos = getPos(i);
        if(_doTiming > 5) _watch->SetTime("initSeedCircle-fitter");
        _circleFitter.addPoint(pos.x(), pos.y(), wP);
        if(_doTiming > 5) _watch->StopTime("initSeedCircle-fitter");
        hit.used = true;
      } else {
        hit.used = false;
      }
    }

    // check if there are enough hits to continue with search
    if (_circleFitter.qn() < _minSeedCircleHits) { outcome = CONTINUE;}
    else { outcome = GOOD; }

    if(_diagLevel > 0) {
      _diagInfo.loopCondition = outcome;
      _hmanager->fillHistograms(&_diagInfo, DIAG::kCircle);
    }
    if(_doTiming > 4) _watch->StopTime(__func__);
  }

  //-----------------------------------------------------------------------------
  // function to initialize phi info relative to helix center in _tcHits
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::initHelixPhi() {
    if(_doTiming > 4) _watch->SetTime(__func__);

    float xC = _circleFitter.x0();
    float yC = _circleFitter.y0();

    for (size_t i = 0; i < _tcHits.size(); i++) {
      if (_tcHits[i].inHelix == true || _tcHits[i].hitIndice == HitType::STOPPINGTARGET || _tcHits[i].used == false) {
        continue;
      }
      // initialize phi data member relative to circle center
      computeHelixPhi(i, xC, yC);
      _tcHits[i].helixPhiCorrection = 0;
    }
    if(_doTiming > 4) _watch->StopTime(__func__);
  }

  //-----------------------------------------------------------------------------
  // function to find seed phi lines
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::findSeedPhiLines(LoopCondition& outcome) {
    if(_doTiming > 4) _watch->SetTime(__func__);

    _seedPhiLines.clear();
    lineSegmentInfo lsInfo;

    const float rC = _circleFitter.radius();

    float lastAddedNegativePhi = 0.f;
    float lastAddedNegativePhiError2 = 0.f;
    float lastAddedNegativeZ = 0.f;
    float lastAddedPositivePhi = 0.f;
    float lastAddedPositivePhiError2 = 0.f;
    float lastAddedPositiveZ = 0.f;
    float maxHitGap = 0.f; // For diagnostics

    // seed point
    lineInfo positiveLine; // for negative/positive slope hypotheses
    lineInfo negativeLine;
    const size_t ntcHits = _tcHits.size();
    const size_t max_index = ntcHits - _minLineSegmentHits; // max + 1
    for (size_t i = 0; i < max_index; i++) {
      auto& tcHit_i = _tcHits[i];
      if (tcHit_i.inHelix == true || tcHit_i.used == false
          || tcHit_i.notOnSegment == false
          || tcHit_i.hitIndice == HitType::STOPPINGTARGET
          || tcHit_i.hitIndice == HitType::CALOCLUSTER) { continue; }

      if(_doTiming > 5) _watch->SetTime("findSeedPhiLines-seed");

      float seedZ = tcHit_i.pos.z();
      float seedPhi = tcHit_i.helixPhi;
      float seedError2 = tcHit_i.helixPhiError2;
      float seedWeight = 1.f / (seedError2);

      // initialize line for both positive and negative slopes
      positiveLine.clear();
      positiveLine.tcHitsIndices.push_back(i);
      positiveLine.helixPhiCorrections.push_back(0);
      positiveLine.zMax = seedZ;
      positiveLine.fitter.addPoint(seedZ, seedPhi, seedWeight);

      negativeLine.clear();
      negativeLine.tcHitsIndices.push_back(i);
      negativeLine.helixPhiCorrections.push_back(0);
      negativeLine.zMax = seedZ;
      negativeLine.fitter.addPoint(seedZ, seedPhi, seedWeight);

      // make clear phi for the most recently added point to fitter
      lastAddedNegativePhi = seedPhi;
      lastAddedNegativePhiError2 = seedError2;
      lastAddedNegativeZ = seedZ;
      lastAddedPositivePhi = seedPhi;
      lastAddedPositivePhiError2 = seedError2;
      lastAddedPositiveZ = seedZ;

      // now loop over test points
      for (size_t j = i + 1; j < ntcHits; j++) {
        auto& tcHit_j = _tcHits[j];
        if (tcHit_j.inHelix == true || tcHit_j.used == false
            || tcHit_j.hitIndice == HitType::STOPPINGTARGET
            || tcHit_j.hitIndice == HitType::CALOCLUSTER) { continue; }

         // Note: looping in decreasing z order, so seed z,last positive z, and last negative z >= test z
        float testZ = tcHit_j.pos.z();
        if (seedZ - testZ > _maxZWindow) { break; }
        float positiveDeltaZ = lastAddedPositiveZ - testZ;
        float negativeDeltaZ = lastAddedNegativeZ - testZ;
        if (positiveDeltaZ > _maxSeedLineGap && negativeDeltaZ > _maxSeedLineGap) { break; }

        float testPhi              = tcHit_j.helixPhi;
        float testPhiError2        = tcHit_j.helixPhiError2;
        float testWeight           = 1.f / (testPhiError2);
        float positiveDiff         = lastAddedPositivePhi - testPhi; // sign dropped in squaring
        float positiveDiffSq       = positiveDiff*positiveDiff; // faster than sqrt(error)
        float positiveDiffErrorSq  = lastAddedPositivePhiError2 + testPhiError2;
        float negativeDiff         = lastAddedNegativePhi - testPhi; // sign dropped in squaring
        float negativeDiffSq       = negativeDiff*negativeDiff; // faster than sqrt(error)
        float negativeDiffErrorSq  = lastAddedNegativePhiError2 + testPhiError2;

        // Test adding the hit to the positive line
        bool add_to_positive = ((lastAddedPositivePhi > testPhi || positiveDiffSq < positiveDiffErrorSq)
                                && positiveDeltaZ <= _maxSeedLineGap && positiveDeltaZ > 1.e-10);
        add_to_positive &= (positiveDiffSq < positiveDiffErrorSq || std::abs((lastAddedPositivePhi-testPhi)/positiveDeltaZ) <= _maxDphiDz
                            || underEstimateSlope(lastAddedPositivePhi,lastAddedPositivePhiError2,testPhi,testPhiError2,positiveDeltaZ) <= _maxDphiDz);
        if(add_to_positive) {
          positiveLine.tcHitsIndices.push_back(j);
          positiveLine.helixPhiCorrections.push_back(0);
          positiveLine.zMin = testZ;
          positiveLine.fitter.addPoint(testZ, testPhi, testWeight);
          lastAddedPositivePhi = testPhi;
          lastAddedPositivePhiError2 = testPhiError2;
          if (_diagLevel > 0 && positiveDeltaZ > maxHitGap) { maxHitGap = positiveDeltaZ; } // for diagnostics, focus on positive lines
        }

        // Test adding the hit to the negative line
        bool add_to_negative = ((testPhi > lastAddedNegativePhi || negativeDiffSq < negativeDiffErrorSq)
                                && negativeDeltaZ <= _maxSeedLineGap && negativeDeltaZ > 1.e-10);
        add_to_negative &= (negativeDiffSq < negativeDiffErrorSq || std::abs((lastAddedNegativePhi-testPhi)/positiveDeltaZ) <= _maxDphiDz
                             || underEstimateSlope(lastAddedNegativePhi,lastAddedNegativePhiError2,testPhi,testPhiError2,negativeDeltaZ) <= _maxDphiDz);
        if(add_to_negative) {
          negativeLine.tcHitsIndices.push_back(j);
          negativeLine.helixPhiCorrections.push_back(0);
          negativeLine.zMin = testZ;
          negativeLine.fitter.addPoint(testZ, testPhi, testWeight);
          lastAddedNegativePhi = testPhi;
          lastAddedNegativePhiError2 = testPhiError2;
        }
      } // end hit j loop

      if (positiveLine.fitter.qn() >= _minLineSegmentHits && positiveLine.fitter.dydx() <= _maxDphiDz) {
        float dphidz = positiveLine.fitter.dydx();
        float p2 = computeHelixMomentum2(rC, dphidz);
        const float dphi = (seedPhi - lastAddedPositivePhi);
        if (dphi*dphi > _segMultiplierSq * seedError2 &&
            p2 > _minHelixMomentumSq &&
            p2 < _maxHelixMomentumSq &&
            positiveLine.fitter.chi2Dof() <= _maxSegmentChi2) {
          // Accept the positive seed phi line
          _seedPhiLines.push_back(positiveLine);
          if (_diagLevel > 0) {
            lsInfo.chi2dof = positiveLine.fitter.chi2Dof();
            lsInfo.maxHitGap = maxHitGap;
            _diagInfo.lineSegmentData.push_back(lsInfo);
          }
          for (size_t n = 0; n < positiveLine.tcHitsIndices.size(); n++) {
            _tcHits[positiveLine.tcHitsIndices[n]].notOnSegment = false;
          }
        }
      }
      if (negativeLine.fitter.qn() >= _minLineSegmentHits && std::abs(negativeLine.fitter.dydx()) <= _maxDphiDz) {
        float dphidz = negativeLine.fitter.dydx();
        float p2 = computeHelixMomentum2(rC, dphidz);
        const float dphi = (seedPhi - lastAddedNegativePhi);
        if (dphi*dphi > _segMultiplierSq * seedError2 &&
            p2 > _minHelixMomentumSq &&
            p2 < _maxHelixMomentumSq &&
            negativeLine.fitter.chi2Dof() <= _maxSegmentChi2) {
          // Accept the negative seed phi line
          _seedPhiLines.push_back(negativeLine);
          for (size_t n = 0; n < negativeLine.tcHitsIndices.size(); n++) {
            _tcHits[negativeLine.tcHitsIndices[n]].notOnSegment = false;
          }
        }
      }
      if(_doTiming > 5) _watch->StopTime("findSeedPhiLines-seed");
    }

    // check if we found any lines
    if (_seedPhiLines.size() == 0) { outcome = CONTINUE; }
    else { outcome = GOOD; }

    if(_diagLevel > 0) {
      _diagInfo.loopCondition = outcome;
      _hmanager->fillHistograms(&_diagInfo, DIAG::kSegments);
    }

    if(_doTiming > 4) _watch->StopTime(__func__);
  }

  //-----------------------------------------------------------------------------
  // function to find seed phi lines
  //-----------------------------------------------------------------------------
  float AgnosticHelixFinder::underEstimateSlope(float phi1, float phi1Err2, float phi2, float phi2Err2, float dz) {

    float slope = 0.0;
    float phi1Err = std::sqrt(phi1Err2);
    float phi2Err = std::sqrt(phi2Err2);

    if (phi2 > phi1) { slope = std::abs(((phi2-phi2Err)-(phi1+phi1Err))/(dz)); }
    else             { slope = std::abs(((phi2+phi2Err)-(phi1-phi1Err))/(dz)); }

    return slope;

  }

  //-----------------------------------------------------------------------------
  // function to resolve 2 Pi Ambiguities for each seed phi line found
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::resolve2PiAmbiguities() {
    if(_doTiming > 4) _watch->SetTime(__func__);

    for (size_t i = 0; i < _seedPhiLines.size(); i++) {
      float lineSlope = _seedPhiLines[i].fitter.dydx();
      float lineIntercept = _seedPhiLines[i].fitter.y0();
      float slopeError = _seedPhiLines[i].fitter.dydxErr();
      float interceptError = _seedPhiLines[i].fitter.y0Err();
      for (size_t j = 0; j < _tcHits.size(); j++) {
        if (_tcHits[j].inHelix == true || _tcHits[j].hitIndice == HitType::STOPPINGTARGET || _tcHits[j].used == false) {
          continue;
        }
        bool alreadyInLine = false;
        float testZ = getPos(j).z();
        float testPhi = _tcHits[j].helixPhi;
        float testPhiError2 = _tcHits[j].helixPhiError2;
        float testWeight = 1.0 / (testPhiError2);
        for (size_t k = 0; k < _seedPhiLines[i].tcHitsIndices.size(); k++) {
          if (_seedPhiLines[i].tcHitsIndices[k] == j) {
            alreadyInLine = true;
            break;
          }
        }
        if (alreadyInLine == true) { continue; }
        // find integer needed for 2pi correction
        float deltaPhi = lineSlope * testZ + lineIntercept - testPhi;
        int n = std::round(deltaPhi / (_twopi));
        testPhi = testPhi + _twopi * n;
        // now see if 2pi corrected phi is consistent with line
        deltaPhi = std::abs(lineSlope * testZ + lineIntercept - testPhi);
        float error2 =
          testZ * testZ * slopeError * slopeError + interceptError * interceptError + testPhiError2;
        if (deltaPhi * deltaPhi < _max2PiAmbigResidual * _max2PiAmbigResidual * error2) {
          _seedPhiLines[i].fitter.addPoint(testZ, testPhi, testWeight);
          _seedPhiLines[i].tcHitsIndices.push_back(j);
          _seedPhiLines[i].helixPhiCorrections.push_back(n);
          if (testZ < _seedPhiLines[i].zMin) { _seedPhiLines[i].zMin = testZ; }
          if (testZ > _seedPhiLines[i].zMax) { _seedPhiLines[i].zMax = testZ; }
        }
      }
    }
    if(_doTiming > 4) _watch->StopTime(__func__);
  }

  //-----------------------------------------------------------------------------
  // function for refining the phi line that was found
  //-----------------------------------------------------------------------------
  bool AgnosticHelixFinder::refinePhiLine(size_t lineIndex) {
    if(_doTiming > 4) _watch->SetTime(__func__);

    float largestResidual2 = 0.0;
    size_t rmIndex = 0;

    for (size_t i = 0; i < _seedPhiLines[lineIndex].tcHitsIndices.size(); i++) {
      size_t tcHitsIndex = _seedPhiLines[lineIndex].tcHitsIndices[i];
      float pointPhi =
        _tcHits[tcHitsIndex].helixPhi + _twopi * _seedPhiLines[lineIndex].helixPhiCorrections[i];
      float pointZ = getPos(tcHitsIndex).z();
      float linePhi =
        _seedPhiLines[lineIndex].fitter.dydx() * pointZ + _seedPhiLines[lineIndex].fitter.y0();
      float deltaPhi = std::abs(linePhi - pointPhi);
      float pointPhiError2 = _tcHits[tcHitsIndex].helixPhiError2;
      float residual2 = deltaPhi * deltaPhi / pointPhiError2;
      if (residual2 > largestResidual2) {
        largestResidual2 = residual2;
        rmIndex = i;
      }
    }

    bool res = false;
    if (largestResidual2 > _maxPhiZResidual * _maxPhiZResidual) {
      size_t tcHitsIndex = _seedPhiLines[lineIndex].tcHitsIndices[rmIndex];
      float z = getPos(tcHitsIndex).z();
      float phi = _tcHits[tcHitsIndex].helixPhi +
        _twopi * _seedPhiLines[lineIndex].helixPhiCorrections[rmIndex];
      float weight = 1.0 / (_tcHits[tcHitsIndex].helixPhiError2);
      _seedPhiLines[lineIndex].fitter.removePoint(z, phi, weight);
      _seedPhiLines[lineIndex].tcHitsIndices.erase(_seedPhiLines[lineIndex].tcHitsIndices.begin() +
                                                   rmIndex);
      _seedPhiLines[lineIndex].helixPhiCorrections.erase(
                                                         _seedPhiLines[lineIndex].helixPhiCorrections.begin() + rmIndex);
      res = true;
    }
    if(_doTiming > 4) _watch->StopTime(__func__);
    return res;
  }

  //-----------------------------------------------------------------------------
  // initialize final circle / line seed prior to hit recovery
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::initFinalSeed(LoopCondition& outcome) {
    if(_doTiming > 4) _watch->SetTime(__func__);

    // make note of number of hits on circle before circle gets updated using the line info
    float seedCircleHits = _circleFitter.qn();

    // first find best phi line by N(hits), and then chi^2/dof for equal N(hits)
    size_t bestLineIndex = 0;
    size_t nHitsMost = 0;
    float chi2dof = 0.0;
    for (size_t i = 0; i < _seedPhiLines.size(); i++) {
      if (_seedPhiLines[i].tcHitsIndices.size() < nHitsMost) { continue; }
      if (_seedPhiLines[i].tcHitsIndices.size() > nHitsMost) {
        nHitsMost = _seedPhiLines[i].tcHitsIndices.size();
        bestLineIndex = i;
        chi2dof = _seedPhiLines[i].fitter.chi2Dof();
        continue;
      }
      if (_seedPhiLines[i].tcHitsIndices.size() == nHitsMost) {
        if (_seedPhiLines[i].fitter.chi2Dof() < chi2dof) {
          bestLineIndex = i;
          chi2dof = _seedPhiLines[i].fitter.chi2Dof();
        }
      }
    }

    // now initialize line/circle fitters and flag hits in tcHits that are used
    _circleFitter.clear();
    _lineFitter.clear();
    _lineFitter = _seedPhiLines[bestLineIndex].fitter;
    for (size_t i = 0; i < _tcHits.size(); i++) {
      if (_tcHits[i].inHelix == true || _tcHits[i].hitIndice == HitType::STOPPINGTARGET) { continue; }
      _tcHits[i].used = false;
      for (size_t j = 0; j < _seedPhiLines[bestLineIndex].tcHitsIndices.size(); j++) {
        if (_seedPhiLines[bestLineIndex].tcHitsIndices[j] == i) {
          _tcHits[i].used = true;
          _tcHits[i].helixPhiCorrection = _seedPhiLines[bestLineIndex].helixPhiCorrections[j];
          float xP = getPos(i).x();
          float yP = getPos(i).y();
          float wP = 1.0 / (_tcHits[i].circleError2);
          _circleFitter.addPoint(xP, yP, wP);
          break;
        }
      }
    }

    // see what ratio of hits from initial seed circle survived the line search
    float nLineHits = _lineFitter.qn();
    float nHitsRatio = seedCircleHits/nLineHits;

    if (_diagLevel > 0) {
      finalLineInfo lInfo;
      lInfo.nHitsRatio = nHitsRatio;
      _diagInfo.lineInfoData.push_back(lInfo);
    }

    // check if we should continue with helix search or not
    if (_circleFitter.qn() < _minFinalSeedHits || nHitsRatio > _maxNHitsRatio) { outcome = CONTINUE;}
    else { outcome = GOOD; }

    if(_diagLevel > 0) {
      _diagInfo.loopCondition = outcome;
      _hmanager->fillHistograms(&_diagInfo, DIAG::kLine);
    }
    if(_doTiming > 4) _watch->StopTime(__func__);
  }

  //-----------------------------------------------------------------------------
  // recover points
  //-----------------------------------------------------------------------------
  bool AgnosticHelixFinder::recoverPoints() {
    if(_doTiming > 4) _watch->SetTime(__func__);

    float xC = _circleFitter.x0();
    float yC = _circleFitter.y0();
    float rC = _circleFitter.radius();

    float lineSlope = _lineFitter.dydx();
    float lineIntercept = _lineFitter.y0();

    float smallestResidual2 = 0.0;
    size_t addIndex = 0;
    bool recoveries = false;

    // first check hits for best recovery
    for (size_t i = 0; i < _tcHits.size(); i++) {
      if (_tcHits[i].used == true || _tcHits[i].inHelix == true || _tcHits[i].hitIndice == HitType::STOPPINGTARGET) {
        continue;
      }
      // compute circle residual
      float circleResidual2 = computeCircleResidual2(i, xC, yC, rC);
      if (circleResidual2 > _maxCircleRecoverSigma * _maxCircleRecoverSigma) {
        continue;
      }
      // compute phi residual
      computeHelixPhi(i, xC, yC);
      float phiP = _tcHits[i].helixPhi;
      float zP = getPos(i).z();
      float phiSigma2 = _tcHits[i].helixPhiError2;
      float phiDistance = lineSlope * zP + lineIntercept - phiP;
      int n = std::round(phiDistance / (_twopi));
      _tcHits[i].helixPhiCorrection = n;
      phiP = phiP + _twopi * n;
      phiDistance = std::abs(lineSlope * zP + lineIntercept - phiP);
      float phiResidual2 = phiDistance * phiDistance / phiSigma2;
      if (phiResidual2 > _maxLineRecoverSigma * _maxLineRecoverSigma) { continue; }
      float residual2 = phiResidual2 + circleResidual2;
      if (residual2 < smallestResidual2 || !recoveries) {
        smallestResidual2 = residual2;
        addIndex = i;
      }
      recoveries = true; // flag that at least one recovery was made
    }

    // add point to circle fitter to update circle parameters
    // recompute circle points weights
    // update circle with new weights
    // recompute phis / phi errors
    // update line
    if (recoveries) {
      _tcHits[addIndex].used = true;
      float w = 1.0 / (_tcHits[addIndex].circleError2);
      _circleFitter.addPoint(getPos(addIndex).x(), getPos(addIndex).y(), w);
      xC = _circleFitter.x0();
      yC = _circleFitter.y0();
      rC = _circleFitter.radius();
      _circleFitter.clear();
      for (size_t i = 0; i < _tcHits.size(); i++) {
        if (_tcHits[i].used == false || _tcHits[i].inHelix == true || _tcHits[i].hitIndice == HitType::STOPPINGTARGET) {
          continue;
        }
        computeCircleError2(i, xC, yC);
        float wP = 1.0 / (_tcHits[i].circleError2);
        _circleFitter.addPoint(getPos(i).x(), getPos(i).y(), wP);
      }
      xC = _circleFitter.x0();
      yC = _circleFitter.y0();
      rC = _circleFitter.radius();
      _lineFitter.clear();
      for (size_t i = 0; i < _tcHits.size(); i++) {
        if (_tcHits[i].used == false || _tcHits[i].inHelix == true || _tcHits[i].hitIndice == HitType::STOPPINGTARGET) {
          continue;
        }
        float zP = getPos(i).z();
        computeHelixPhi(i, xC, yC);
        float phiWeight = 1.0 / (_tcHits[i].helixPhiError2);
        float phiP = _tcHits[i].helixPhi;
        float deltaPhi = lineSlope * zP + lineIntercept - phiP;
        _tcHits[i].helixPhiCorrection = std::round(deltaPhi / (_twopi));
        phiP = phiP + _tcHits[i].helixPhiCorrection * _twopi;
        _lineFitter.addPoint(zP, phiP, phiWeight);
      }
    }
    if(_diagLevel > 0) _hmanager->fillHistograms(&_diagInfo, DIAG::kRecover);
    if(_doTiming > 4) _watch->StopTime(__func__);
    return recoveries;
  }

  //-----------------------------------------------------------------------------
  // check if helix is savable based on fcl defined cuts
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::checkHelixViability(LoopCondition& outcome) {
    if(_doTiming > 4) _watch->SetTime(__func__);

    // first check if line has good enough fit
    if (_lineFitter.chi2Dof() > _chi2LineSaveThresh) {
      outcome = CONTINUE;
    } else {
      // if line fit was good enough make sure momenta are within allowable range
      float radius = _circleFitter.radius();
      float slope = _lineFitter.dydx();
      float mom2 = computeHelixMomentum2(radius, slope);
      float pt = computeHelixPerpMomentum(radius);
      if (mom2 < _minHelixMomentumSq ||
          mom2 > _maxHelixMomentumSq || pt < _minHelixPerpMomentum ||
          pt > _maxHelixPerpMomentum) { outcome = CONTINUE; }
      else { outcome = GOOD; }
    }

    if(_diagLevel > 0) {
      _diagInfo.loopCondition = outcome;
      _hmanager->fillHistograms(&_diagInfo, DIAG::kViability);
    }
    if(_doTiming > 4) _watch->StopTime(__func__);
  }

  //-----------------------------------------------------------------------------
  // function to save helix
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::saveHelix(size_t tc, HelixSeedCollection& HSColl) {
    HelixSeed hseed;
    hseed._t0 = _tcColl->at(tc)._t0;
    hseed._timeCluster = art::Ptr<mu2e::TimeCluster>(_tcCollH, tc);
    hseed._hhits.setParent(_chColl->parent());

    // flag hits used in helix, and push to combo hit collection in helix seed
    // also add points to linear fitter to get t0
    ::LsqSums2 fitter;
    for (size_t i = 0; i < _tcHits.size(); i++) {
      if (_tcHits[i].inHelix == true || _tcHits[i].used == false
          || _tcHits[i].hitIndice == HitType::STOPPINGTARGET
          || _tcHits[i].hitIndice == HitType::CALOCLUSTER) { continue; }
      _tcHits[i].inHelix = true;
      const ComboHit& hit = *(_tcHits[i].hit);
      fitter.addPoint(hit.pos().z(), hit.correctedTime(), 1 / (hit.timeRes() * hit.timeRes()));
      ComboHit hhit(hit);
      hhit._hphi = _tcHits[i].helixPhi + _tcHits[i].helixPhiCorrection * _twopi;
      hseed._hhits.push_back(hhit);
    }
    float eDepAvg = hseed._hhits.eDepAvg();

    hseed._t0 = TrkT0(fitter.y0(), fitter.y0Err());

    // set the geom helix parameters
    float x0 = _circleFitter.x0();
    float y0 = _circleFitter.y0();
    CLHEP::Hep3Vector center(x0, y0, 0);
    float dfdz = _lineFitter.dydx();

    hseed._helix._fz0 = std::fmod(_lineFitter.y0(), _twopi);
    if     (hseed._helix._fz0 >  M_PI) { hseed._helix._fz0 = hseed._helix._fz0 - _twopi; }
    else if(hseed._helix._fz0 < -M_PI) { hseed._helix._fz0 = hseed._helix._fz0 + _twopi; }
    hseed._helix._rcent = center.perp();
    hseed._helix._fcent = center.phi();
    hseed._helix._radius = _circleFitter.radius();
    hseed._helix._lambda = 1. / dfdz;
    hseed._helix._helicity = hseed._helix._lambda > 0 ? Helicity::poshel : Helicity::neghel;

    // include also the values of the chi2
    hseed._helix._chi2dXY = _circleFitter.chi2DofCircle();
    hseed._helix._chi2dZPhi = _lineFitter.chi2Dof();

    hseed._status.merge(TrkFitFlag::helixOK);
    hseed._status.merge(TrkFitFlag::APRHelix);

    if(_diagLevel > 0) { // diagnostics for all helix seeds, independent of success
      _diagInfo.hseed = &hseed;
      _hmanager->fillHistograms(&_diagInfo, DIAG::kHelix);
    }

    if (eDepAvg > _maxEDepAvg) return;

    // compute direction of propagation and make save decision
    HelixTool helTool(&hseed,_tracker);
    float tzSlope = 0.0;
    float tzSlopeErr = 0.0;
    float tzSlopeChi2 = 0.0;
    helTool.dirOfProp(tzSlope, tzSlopeErr, tzSlopeChi2);
    HelixRecoDir helDir(tzSlope, tzSlopeErr, tzSlopeChi2);
    hseed._recoDir = helDir;
    hseed._propDir = helDir.predictDirection(_tzSlopeSigThresh);
    if (!validHelixDirection(hseed._propDir)) return;

    // push back the helix seed to the helix seed collection
    HSColl.emplace_back(hseed);

    if(_diagLevel > 0) { // diagnostics for accepted helix seeds
      _diagInfo.helixSeedData.push_back(hsInfo(_bz0, hseed)); // list of helices per time cluster
      _hmanager->fillHistograms(&_diagInfo, DIAG::kFinal);
    }

  }

  //-----------------------------------------------------------------------------
  // check if the propagation direction is among the desired save directions
  //-----------------------------------------------------------------------------
  bool AgnosticHelixFinder::validHelixDirection(TrkFitDirection::FitDirection direction) {

    for (size_t i=0; i<_validHelixDirections.size(); i++) {
      if (_validHelixDirections.at(i) == direction) return true;
    }

    return false;

  }

} // namespace mu2e

using mu2e::AgnosticHelixFinder;
DEFINE_ART_MODULE(AgnosticHelixFinder)
