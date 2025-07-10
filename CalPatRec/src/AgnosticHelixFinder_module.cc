///////////////////////////////////////////////////////////////////////////////
// AgnosticHelixFinder
// M. Stortini & E. Martinez
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

    struct cHit {
      int     hitIndice = 0; // index of point in _chColl
      float   circleError2 = 1.0;
      float   helixPhi = 0.0;
      float   helixPhiError2 = 0.0;
      int     helixPhiCorrection = 0;
      bool    inHelix = false;
      bool    used = false; // whether or not hit is used in fits
      bool    isolated = false;
      bool    averagedOut = false;
      bool    highEDep = false; // avoid hits in tripleting that are most likely protons
      bool    notOnLine = true;
      bool    uselessTripletSeed = false;
      bool    notOnSegment = true;
    };

    struct tripletPoint {
      XYZVectorF   pos;
      int          hitIndice;
    };

    struct triplet {
      tripletPoint i;
      tripletPoint j;
      tripletPoint k;
    };

    struct lineInfo {
      float                 zMin;
      float                 zMax;
      ::LsqSums2            fitter;
      std::vector<size_t>   tcHitsIndices;
      std::vector<int>      helixPhiCorrections; // integer for 2pi ambiguity
    };

    enum LoopCondition {
      CONTINUE,
      BREAK,
      GOOD
    };

    enum HitType {
      CALOCLUSTER = -1,
      STOPPINGTARGET = -2
    };

  private:
    //-----------------------------------------------------------------------------
    // tool on or off
    //-----------------------------------------------------------------------------
    int _diagLevel;

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
    int      _isoMinHitsNear;
    bool     _doAverageFlag;
    bool     _doEDepFlag;
    float    _eDepFlagThresh;
    float    _minDistCut;
    float    _minTripletSeedZ;
    float    _minTripletDz;
    float    _maxTripletDz;
    float    _minTripletDist;
    float    _minTripletArea;
    float    _maxSeedCircleResidual;
    int      _minSeedCircleHits;
    float    _maxDphiDz;
    float    _maxSeedLineGap;
    float    _maxZWindow;
    int      _minLineSegmentHits;
    float    _segMultiplier;
    float    _maxSegmentChi2;
    float    _max2PiAmbigResidual;
    float    _maxPhiZResidual;
    int      _minFinalSeedHits;
    float    _maxNHitsRatio;
    float    _maxCircleRecoverSigma;
    float    _maxLineRecoverSigma;
    float    _caloClusterSigma;
    int      _minNHelixStrawHits;
    int      _minNHelixComboHits;
    float    _minHelixPerpMomentum;
    float    _maxHelixPerpMomentum;
    float    _minHelixMomentum;
    float    _maxHelixMomentum;
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
    art::Event*                          _event;

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
    static constexpr float mmTconversion = CLHEP::c_light/1000.0;

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
    XYZVectorF   getPos                    (size_t& tcHitsIndex);
    void         computeCircleError2       (size_t& tcHitsIndex, float& xC, float& yC);
    float        computeCircleResidual2    (size_t& tcHitsIndex, float& xC, float& yC, float& rC);
    void         computeHelixPhi           (size_t& tcHitsIndex, float& xC, float& yC);
    float        computeHelixMomentum2     (float& radius, float& dphidz);
    float        computeHelixPerpMomentum  (float& radius);
    void         tcHitsFill                (size_t tc);
    void         setFlags                  ();
    void         resetFlags                ();
    void         findHelix                 (size_t tc, HelixSeedCollection& HSColl, bool& findAnotherHelix);
    bool         passesFlags               (size_t& tcHitsIndex);
    void         setTripletI               (size_t& tcHitsIndex, triplet& trip, LoopCondition& outcome);
    void         setTripletJ               (size_t& tcHitsIndex, triplet& trip, LoopCondition& outcome);
    void         setTripletK               (size_t& tcHitsIndex, triplet& trip, LoopCondition& outcome);
    void         initTriplet               (triplet& trip, LoopCondition& outcome);
    void         initSeedCircle            (LoopCondition& outcome);
    void         initHelixPhi              ();
    void         findSeedPhiLines          (LoopCondition& outcome);
    float        underEstimateSlope        (float& phi1, float& phi1Err2, float& phi2, float& phi2Err2, float& dz);
    void         resolve2PiAmbiguities     ();
    void         refinePhiLine             (size_t lineIndex, bool& removals);
    void         initFinalSeed             (LoopCondition& outcome);
    void         recoverPoints             (bool& recoveries);
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

    // convert the helix direction names into enums
    for(auto helix_dir : config().validHelixDirections()) {
      _validHelixDirections.push_back(TrkFitDirection::fitDirectionFromName(helix_dir));
    }
    consumes<ComboHitCollection>     (_chLabel);
    consumes<TimeClusterCollection>  (_tcLabel);
    consumes<CaloClusterCollection>  (_ccLabel);
    produces<HelixSeedCollection>    ();

    if (_diagLevel == 1) _hmanager = art::make_tool<ModuleHistToolBase>(config().diagPlugin, "diagPlugin");
    else _hmanager = std::make_unique<ModuleHistToolBase>();

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
    if (_diagLevel == 1) {
      art::ServiceHandle<art::TFileService> tfs;
      _hmanager->bookHistograms(tfs);
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
  }

  //-----------------------------------------------------------------------------
  // find input things
  //-----------------------------------------------------------------------------
  bool AgnosticHelixFinder::findData(const art::Event& evt) {

    auto chCollH = evt.getValidHandle<ComboHitCollection>(_chLabel);
    if (chCollH.product() != 0) {
      _chColl = chCollH.product();
    } else { _chColl = 0; }

    auto _tcCollH = evt.getValidHandle<TimeClusterCollection>(_tcLabel);
    if (_tcCollH.product() != 0) {
      _tcColl = _tcCollH.product();
    } else { _tcColl = 0; }

    if (evt.getByLabel(_ccLabel, _ccHandle)) {
      _ccColl = _ccHandle.product();
    } else { _ccColl = NULL; }

    return (_tcColl != 0);
  }

  //-----------------------------------------------------------------------------
  // event entry point
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::produce(art::Event& event) {

    // set the event and prepare hsColl
    _event = &event;
    std::unique_ptr<HelixSeedCollection> hsColl(new HelixSeedCollection);

    // get the necessary data to do helix search
    bool dataExists = findData(event);

    // prepare diagnostic tool data members
    if (_diagLevel == 1) {
      _diagInfo.nHelices = 0;
      _diagInfo.timeClusterData.clear();
      _diagInfo.helixSeedData.clear();
      _diagInfo.lineSegmentData.clear();
    }

    // flag whether event is intense or not
    if ((int)_tcColl->size() > _intenseEventThresh) {
      _intenseEvent = true;
    } else { _intenseEvent = false; }

    // do helix search
    bool continueSearch = true;

    // search for helix
    if (dataExists) {
      for (size_t i = 0; i < _tcColl->size(); i++) {
        // check to see if cluster is a busy one
        _intenseCluster = false;
        // if cluster is busy or event is intense then only process time cluster if it has calo cluster
        if ((int)_tcColl->at(i).nhits() > _intenseClusterThresh) { _intenseCluster = true; }
        if (_intenseEvent == true || _intenseCluster == true) {
          if (_tcColl->at(i).hasCaloCluster() == false) { continue; }
        }
        tcInfo timeClusterInfo;
        int nHelicesInitial = _diagInfo.nHelices;
        _tcHits.clear();
        tcHitsFill(i);
        continueSearch = true;
        while (continueSearch == true) {
          findHelix(i, *hsColl, continueSearch);
          if (_findMultipleHelices == false) { continueSearch = false; }
        }
        if (_diagLevel == 1) {
          timeClusterInfo.nHelices = _diagInfo.nHelices - nHelicesInitial;
          timeClusterInfo.nComboHits = _tcColl->at(i).nhits();
          timeClusterInfo.nStrawHits = _tcColl->at(i).nStrawHits();
          _diagInfo.timeClusterData.push_back(timeClusterInfo);
        }
      }
    }

    // put helix seed collection into the event record
    event.put(std::move(hsColl));

    // fill necessary data members for diagnostic tool
    if (_diagLevel == 1) {
      _diagInfo.nTimeClusters = _tcColl->size();
      _hmanager->fillHistograms(&_diagInfo);
    }
  }

  //-----------------------------------------------------------------------------
  // endJob
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::endJob() {}

  //-----------------------------------------------------------------------------
  // get pointer to position of hit given some index in _tcHits
  //-----------------------------------------------------------------------------
  XYZVectorF AgnosticHelixFinder::getPos(size_t& tcHitsIndex) {

    int hitIndice = _tcHits[tcHitsIndex].hitIndice;

    if (hitIndice >= 0) { return _chColl->at(hitIndice).pos(); }
    if (hitIndice == HitType::STOPPINGTARGET) { return _stopTargPos; }
    if (hitIndice == HitType::CALOCLUSTER) { return _caloPos; }

    XYZVectorF p(0.0,0.0,0.0);
    return p;

  }

  //-----------------------------------------------------------------------------
  // updates circleError of hit given some circle parameters
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::computeCircleError2(size_t& tcHitsIndex, float& xC, float& yC) {

    int hitIndice = _tcHits[tcHitsIndex].hitIndice;

    if (hitIndice == HitType::CALOCLUSTER) {
      _tcHits[tcHitsIndex].circleError2 = _caloClusterSigma * _caloClusterSigma;
    } else {
      float transVar = _chColl->at(hitIndice).transVar();
      float x = getPos(tcHitsIndex).x();
      float y = getPos(tcHitsIndex).y();
      float dx = x - xC;
      float dy = y - yC;
      float dxn = dx * _chColl->at(hitIndice).vDir().x() + dy * _chColl->at(hitIndice).vDir().y();
      float costh2 = dxn * dxn / (dx * dx + dy * dy);
      float sinth2 = 1 - costh2;
      _tcHits[tcHitsIndex].circleError2 =
        _chColl->at(hitIndice).wireVar()*sinth2 + transVar*costh2;
    }
  }

  //-----------------------------------------------------------------------------
  // compute residual between point an circle
  //-----------------------------------------------------------------------------
  float AgnosticHelixFinder::computeCircleResidual2(size_t& tcHitsIndex, float& xC, float& yC, float& rC) {

    float xP = getPos(tcHitsIndex).x();
    float yP = getPos(tcHitsIndex).y();
    float deltaDistance = std::abs(rC - std::sqrt((xP - xC) * (xP - xC) + (yP - yC) * (yP - yC)));
    float circleSigma2 = _tcHits[tcHitsIndex].circleError2;

    return deltaDistance * deltaDistance / circleSigma2;
  }

  //-----------------------------------------------------------------------------
  // compute phi relative to helix center, and set helixPhiError
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::computeHelixPhi(size_t& tcHitsIndex, float& xC, float& yC) {

    int hitIndice = _tcHits[tcHitsIndex].hitIndice;

    float X = getPos(tcHitsIndex).x() - xC;
    float Y = getPos(tcHitsIndex).y() - yC;
    _tcHits[tcHitsIndex].helixPhi = polyAtan2(Y, X);
    if (_tcHits[tcHitsIndex].helixPhi < 0) {
      _tcHits[tcHitsIndex].helixPhi = _tcHits[tcHitsIndex].helixPhi + 2 * M_PI;
    }
    // find phi error and initialize it
    // find phi error by projecting errors onto vector tangent to circle
    float deltaS2(0);
    if (hitIndice == HitType::CALOCLUSTER) {
      deltaS2 = _caloClusterSigma * _caloClusterSigma;
    } else {
      float tanVecX = Y / std::sqrt(X * X + Y * Y);
      float tanVecY = -X / std::sqrt(X * X + Y * Y);
      float wireErr = _chColl->at(hitIndice).wireRes();
      float wireVecX = _chColl->at(hitIndice).uDir().x();
      float wireVecY = _chColl->at(hitIndice).uDir().y();
      float projWireErr = wireErr * (wireVecX * tanVecX + wireVecY * tanVecY);
      float transErr = _chColl->at(hitIndice).transRes();
      float transVecX = _chColl->at(hitIndice).uDir().y();
      float transVecY = -_chColl->at(hitIndice).uDir().x();
      float projTransErr = transErr * (transVecX * tanVecX + transVecY * tanVecY);
      deltaS2 = projWireErr * projWireErr + projTransErr * projTransErr;
    }
    _tcHits[tcHitsIndex].helixPhiError2 = deltaS2 / (X * X + Y * Y);
  }

  //-----------------------------------------------------------------------------
  // function to compute helix momentum^2 given some circle radius and line slope
  //-----------------------------------------------------------------------------
  float AgnosticHelixFinder::computeHelixMomentum2(float& radius, float& dphidz) {

    float lambda = 1.0 / dphidz;

    return (radius * radius + lambda * lambda) * _bz0 * _bz0 * mmTconversion * mmTconversion;
  }

  //-----------------------------------------------------------------------------
  // function to compute helix transverse momentum given circle radius and b-field
  //-----------------------------------------------------------------------------
  float AgnosticHelixFinder::computeHelixPerpMomentum(float& radius) {

    return radius * _bz0 * mmTconversion;
  }

  //-----------------------------------------------------------------------------
  // fill vector with hits to search for helix
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::tcHitsFill(size_t tc) {

    size_t sortStartIndex = 0;

    // push back stopping target if it is to be used
    if (_useStoppingTarget == true) {
      cHit hit;
      hit.hitIndice = HitType::STOPPINGTARGET;
      _tcHits.push_back(hit);
      sortStartIndex++;
    }

    // push back calo cluster if it exists in time cluster
    if (_tcColl->at(tc).hasCaloCluster()) {
      const art::Ptr<CaloCluster> cl = _tcColl->at(tc).caloCluster();
      const CaloCluster*cluster = cl.get();
      if (cluster != nullptr){
        cHit hit;
        hit.hitIndice = HitType::CALOCLUSTER;
        CLHEP::Hep3Vector gpos = _calorimeter->geomUtil().diskToMu2e(cl->diskID(), cl->cog3Vector());
        CLHEP::Hep3Vector tpos = _calorimeter->geomUtil().mu2eToTracker(gpos);
        double offset = _calorimeter->caloInfo().getDouble("diskCaseZLength");
        offset += _calorimeter->caloInfo().getDouble("BPPipeZOffset");
        offset += _calorimeter->caloInfo().getDouble("BPHoleZLength");
        offset += _calorimeter->caloInfo().getDouble("FEEZLength");
        offset /= 2.0;
        _caloPos.SetCoordinates(tpos.x(), tpos.y(), tpos.z() - offset);
        _tcHits.push_back(hit);
        sortStartIndex++;
      }
    }

    // fill hits from time cluster
    for (size_t i = 0; i < _tcColl->at(tc)._strawHitIdxs.size(); i++) {
      cHit hit;
      hit.hitIndice = _tcColl->at(tc)._strawHitIdxs[i];
      _tcHits.push_back(hit);
    }

    // order from largest z to smallest z (skip over stopping target and calo cluster since they
    // aren't in _chColl)
    std::sort(_tcHits.begin() + sortStartIndex, _tcHits.end(), [&](const cHit& a, const cHit& b) {
      return _chColl->at(a.hitIndice).pos().z() > _chColl->at(b.hitIndice).pos().z();
    });
  }

  //-----------------------------------------------------------------------------
  // logic for setting certain flags on hits
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::setFlags() {

    // do isolation, average, and eDepFlag flagging
    if (_doIsolationFlag == true || _doAverageFlag == true || _doEDepFlag == true) {
      for (size_t i = 0; i < _tcHits.size(); i++) {
        if (_tcHits[i].inHelix == true || _tcHits[i].hitIndice < 0) { continue; }
        int hitIndice = _tcHits[i].hitIndice;
        float hitEnergy = _chColl->at(hitIndice).energyDep();
        if (_doEDepFlag == true && hitEnergy > _eDepFlagThresh) { _tcHits[i].highEDep = true; }
        if (_doIsolationFlag == true || _doAverageFlag == true) {
          int nHitsNear = 0;
          XYZVectorF seedPos = getPos(i);
          for (size_t j = 0; j < _tcHits.size(); j++) {
            if (_tcHits[j].inHelix == true || _tcHits[j].hitIndice < 0) { continue; }
            if (j == i) { continue; }
            XYZVectorF testPos = getPos(j);
            // do isolation flagging
            if (_doIsolationFlag == true) {
              if ((seedPos-testPos).Perp2() < _isoRad * _isoRad) { nHitsNear++; }
              if (nHitsNear < _isoMinHitsNear) { _tcHits[i].isolated = true; }
              else { _tcHits[i].isolated = false; }
            }
            // do averaging out
            if (_doAverageFlag == true) {
              if (_tcHits[i].averagedOut == true) { continue; }
              if (_tcHits[j].averagedOut == true) { continue; }
              if ((seedPos-testPos).Perp2() <= _minDistCut * _minDistCut) { _tcHits[j].averagedOut = true; }
           }
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
  void AgnosticHelixFinder::findHelix(size_t tc, HelixSeedCollection& HSColl, bool& findAnotherHelix) {

    // clear line and circle fitters
    _circleFitter.clear();
    _lineFitter.clear();

    // by default we set findAnotherHelix to false, if we find a helix then this will be set to true
    findAnotherHelix = false;

    // set flags before starting search so that we know what hits to use for tripletting
    if (_doIsolationFlag == true || _doAverageFlag == true) { setFlags(); }

    // now we loop over triplets
    bool foundPhiZRemoval = true;
    bool pointRecovered = true;
    for (size_t i = 0; i < _tcHits.size() - 2; i++) {
      bool uselessSeed = true;
      triplet tripletInfo;
      LoopCondition loopCondition;
      setTripletI(i, tripletInfo, loopCondition);
      if (loopCondition == CONTINUE) { continue; }
      if (loopCondition == BREAK) { break; }
      for (size_t j = i + 1; j < _tcHits.size() - 1; j++) {
        setTripletJ(j, tripletInfo, loopCondition);
        if (loopCondition == CONTINUE) { continue; }
        if (loopCondition == BREAK) { break; }
        for (size_t k = j + 1; k < _tcHits.size(); k++) {
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
            foundPhiZRemoval = true;
            while (foundPhiZRemoval == true) {
              refinePhiLine(ii, foundPhiZRemoval);
            }
          }
          initFinalSeed(loopCondition);
          if (loopCondition == CONTINUE) { continue; }
          pointRecovered = true;
          while (pointRecovered == true) { recoverPoints(pointRecovered); }
          checkHelixViability(loopCondition);
          if (loopCondition == CONTINUE) { continue; }
          // before saving helix we make sure it has enough hits
          int nStrawHitsInHelix = 0;
          int nComboHitsInHelix = 0;
          int nStrawHitsInTimeCluster = 0;
          int nComboHitsInTimeCluster = 0;
          // compute number of usable hits in time cluster, and number of hits in candidate helix
          for (size_t q = 0; q < _tcHits.size(); q++) {
            if (_tcHits[q].inHelix == true || _tcHits[q].hitIndice < 0) { continue; }
            int hitIndice = _tcHits[q].hitIndice;
            nStrawHitsInTimeCluster = nStrawHitsInTimeCluster + _chColl->at(hitIndice).nStrawHits();
            nComboHitsInTimeCluster = nComboHitsInTimeCluster + 1;
            if (_tcHits[q].used == false) { continue; }
            nStrawHitsInHelix = nStrawHitsInHelix + _chColl->at(hitIndice).nStrawHits();
            nComboHitsInHelix = nComboHitsInHelix + 1;
          }
          if (nStrawHitsInHelix >= _minNHelixStrawHits && nComboHitsInHelix >= _minNHelixComboHits) {
            saveHelix(tc, HSColl);
            if (_diagLevel == 1) { _diagInfo.nHelices++; }
            // we only want to search for another helix if we have enough remaining hits after saving helix
            int remainingStrawHits = nStrawHitsInTimeCluster - nStrawHitsInHelix;
            int remainingComboHits = nComboHitsInTimeCluster - nComboHitsInHelix;
            if (remainingStrawHits < _minNHelixStrawHits || remainingComboHits < _minNHelixComboHits) {
              return;
            } else {
              findAnotherHelix = true;
              if (_doAverageFlag == true) { resetFlags(); }
              return;
            }
          }
        }
      }
      _tcHits[i].uselessTripletSeed = uselessSeed;
    }
  }

  //-----------------------------------------------------------------------------
  // check flags to see if point is good for triplet-ing with
  //-----------------------------------------------------------------------------
  bool AgnosticHelixFinder::passesFlags(size_t& tcHitsIndex) {

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
  void AgnosticHelixFinder::setTripletI(size_t& tcHitsIndex, triplet& trip, LoopCondition& outcome) {

    // set info for point
    trip.i.pos = getPos(tcHitsIndex);
    trip.i.hitIndice = _tcHits[tcHitsIndex].hitIndice;

    // check if point should break for loop
    if (trip.i.pos.z() < _minTripletSeedZ) {
      outcome = BREAK;
      return;
    }
    if (_intenseEvent == true || _intenseCluster == true) {
      if (trip.i.hitIndice >= 0) {
        outcome = BREAK;
        return;
      }
    }

    // check if point should be continued on
    if (!passesFlags(tcHitsIndex)) {
      outcome = CONTINUE;
      return;
    }

    outcome = GOOD;
  }

  //-----------------------------------------------------------------------------
  // set jth point of triplet, return outcome value which can be used to direct for loops
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::setTripletJ(size_t& tcHitsIndex, triplet& trip, LoopCondition& outcome) {

    // set info for point
    trip.j.pos = getPos(tcHitsIndex);
    trip.j.hitIndice = _tcHits[tcHitsIndex].hitIndice;
    float dz12 = trip.i.pos.z() - trip.j.pos.z();

    // check if point should break for loop
    if (trip.i.hitIndice >= 0 && dz12 > _maxTripletDz) {
      outcome = BREAK;
      return;
    }
    if (_intenseEvent == true || _intenseCluster == true) {
      if (trip.i.hitIndice == HitType::STOPPINGTARGET && trip.j.hitIndice >= 0) {
        outcome = BREAK;
        return;
      }
    }

    // check if point should be continued on
    if (!passesFlags(tcHitsIndex) || dz12 < _minTripletDz ||
        (trip.i.pos-trip.j.pos).Perp2() < _minTripletDist * _minTripletDist) {
      outcome = CONTINUE;
      return;
    }

    outcome = GOOD;
  }

  //-----------------------------------------------------------------------------
  // set kth point of triplet, return outcome value which can be used to direct for loops
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::setTripletK(size_t& tcHitsIndex, triplet& trip, LoopCondition& outcome) {

    // set info for point
    trip.k.pos = getPos(tcHitsIndex);
    trip.k.hitIndice = _tcHits[tcHitsIndex].hitIndice;
    float dz23 = trip.j.pos.z() - trip.k.pos.z();

    // check if point should break for loop
    if (trip.j.hitIndice >= 0 && dz23 > _maxTripletDz) {
      outcome = BREAK;
      return;
    }

    // check if point should be continued on
    if (!passesFlags(tcHitsIndex) || dz23 < _minTripletDz ||
        (trip.i.pos-trip.k.pos).Perp2() < _minTripletDist * _minTripletDist ||
        (trip.j.pos-trip.k.pos).Perp2() < _minTripletDist * _minTripletDist) {
      outcome = CONTINUE;
      return;
    }

    outcome = GOOD;
  }

  //-----------------------------------------------------------------------------
  // finding circle from triplet
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::initTriplet(triplet& trip, LoopCondition& outcome) {

    _circleFitter.addPoint(trip.i.pos.x(), trip.i.pos.y());
    _circleFitter.addPoint(trip.j.pos.x(), trip.j.pos.y());
    _circleFitter.addPoint(trip.k.pos.x(), trip.k.pos.y());

    // check if circle is valid for search or if we should continue to next triplet
    float radius = _circleFitter.radius();
    float pt = computeHelixPerpMomentum(radius);
    if (pt < _minHelixPerpMomentum || pt > _maxHelixPerpMomentum) { outcome = CONTINUE; }
    else { outcome = GOOD; }
  }

  //-----------------------------------------------------------------------------
  // start with initial seed circle
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::initSeedCircle(LoopCondition& outcome) {

    // get triplet circle parameters then clear fitter
    float xC = _circleFitter.x0();
    float yC = _circleFitter.y0();
    float rC = _circleFitter.radius();
    _circleFitter.clear();

    // project error bars onto the triplet circle found and add to fitter those within defined max
    // residual
    for (size_t i = 0; i < _tcHits.size(); i++) {
      if (_tcHits[i].inHelix == true) { continue; }
      if (_tcHits[i].hitIndice == HitType::STOPPINGTARGET) {
        _tcHits[i].used = false;
        continue;
      }
      computeCircleError2(i, xC, yC);
      if (computeCircleResidual2(i, xC, yC, rC) < _maxSeedCircleResidual * _maxSeedCircleResidual) {
        float wP = 1.0 / (_tcHits[i].circleError2);
        _circleFitter.addPoint(getPos(i).x(), getPos(i).y(), wP);
        _tcHits[i].used = true;
      } else {
        _tcHits[i].used = false;
      }
    }

    // check if there are enough hits to continue with search
    if (_circleFitter.qn() < _minSeedCircleHits) { outcome = CONTINUE;}
    else { outcome = GOOD; }
  }

  //-----------------------------------------------------------------------------
  // function to initialize phi info relative to helix center in _tcHits
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::initHelixPhi() {

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
  }

  //-----------------------------------------------------------------------------
  // function to find seed phi lines
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::findSeedPhiLines(LoopCondition& outcome) {

    _seedPhiLines.clear();
    lineSegmentInfo lsInfo;

    float rC = _circleFitter.radius();

    float lastAddedNegativePhi = 0.0;
    float lastAddedNegativePhiError2 = 0.0;
    float lastAddedNegativeZ = 0.0;
    float lastAddedPositivePhi = 0.0;
    float lastAddedPositivePhiError2 = 0.0;
    float lastAddedPositiveZ = 0.0;
    float maxHitGap = 0.0;

    // seed point
    for (size_t i = 0; i < _tcHits.size() - _minLineSegmentHits; i++) {
      if (_tcHits[i].inHelix == true || _tcHits[i].hitIndice < 0) { continue; }
      if (_tcHits[i].used == false) { continue; }
      if (_tcHits[i].notOnSegment == false) { continue; }
      float seedZ = getPos(i).z();
      float seedPhi = _tcHits[i].helixPhi;
      float seedError2 = _tcHits[i].helixPhiError2;
      float seedWeight = 1.0 / (seedError2);
      // initialize line for both positive and negative slopes
      lineInfo _positiveLine;
      _positiveLine.tcHitsIndices.push_back(i);
      _positiveLine.helixPhiCorrections.push_back(0);
      _positiveLine.zMax = seedZ;
      _positiveLine.fitter.addPoint(seedZ, seedPhi, seedWeight);
      lineInfo _negativeLine;
      _negativeLine.tcHitsIndices.push_back(i);
      _negativeLine.helixPhiCorrections.push_back(0);
      _negativeLine.zMax = seedZ;
      _negativeLine.fitter.addPoint(seedZ, seedPhi, seedWeight);
      // make clear phi for the most recently added point to fitter
      lastAddedNegativePhi = seedPhi;
      lastAddedNegativePhiError2 = seedError2;
      lastAddedNegativeZ = seedZ;
      lastAddedPositivePhi = seedPhi;
      lastAddedPositivePhiError2 = seedError2;
      lastAddedPositiveZ = seedZ;
      // now loop over test point
      for (size_t j = i + 1; j < _tcHits.size(); j++) {
        if (_tcHits[j].inHelix == true || _tcHits[j].hitIndice < 0) { continue; }
        if (_tcHits[j].used == false) { continue; }
        float testZ = getPos(j).z();
        if (testZ == seedZ) { continue; }
        if (std::abs(seedZ - testZ) > _maxZWindow) { break; }
        float positiveDeltaZ = lastAddedPositiveZ - testZ;
        float negativeDeltaZ = lastAddedNegativeZ - testZ;
        if (positiveDeltaZ > _maxSeedLineGap && negativeDeltaZ > _maxSeedLineGap) { break; }
        float testPhi = _tcHits[j].helixPhi;
        float testError2 = _tcHits[j].helixPhiError2;
        float testWeight = 1.0 / (testError2);
        float positiveDiff = std::abs(lastAddedPositivePhi - testPhi);
        float positiveDiffError = std::sqrt(lastAddedPositivePhiError2 + testError2);
        float negativeDiff = std::abs(lastAddedNegativePhi - testPhi);
        float negativeDiffError = std::sqrt(lastAddedNegativePhiError2 + testError2);
        if ((lastAddedPositivePhi > testPhi || positiveDiff < positiveDiffError) && positiveDeltaZ <= _maxSeedLineGap) {
          if (positiveDiff >= positiveDiffError && std::abs((lastAddedPositivePhi-testPhi)/positiveDeltaZ) > _maxDphiDz) {
            if (underEstimateSlope(lastAddedPositivePhi,lastAddedPositivePhiError2,testPhi,testError2,positiveDeltaZ) > _maxDphiDz) {
              continue;
            }
          }
          _positiveLine.tcHitsIndices.push_back(j);
          _positiveLine.helixPhiCorrections.push_back(0);
          _positiveLine.zMin = testZ;
          _positiveLine.fitter.addPoint(testZ, testPhi, testWeight);
          lastAddedPositivePhi = testPhi;
          lastAddedPositivePhiError2 = testError2;
          if (positiveDeltaZ > maxHitGap) { maxHitGap = positiveDeltaZ; }
        }
        if ((testPhi > lastAddedNegativePhi || negativeDiff < negativeDiffError) && negativeDeltaZ <= _maxSeedLineGap) {
          if (negativeDiff >= negativeDiffError && std::abs((lastAddedNegativePhi-testPhi)/positiveDeltaZ) > _maxDphiDz) {
            if (underEstimateSlope(lastAddedNegativePhi,lastAddedNegativePhiError2,testPhi,testError2,negativeDeltaZ) > _maxDphiDz) {
              continue;
            }
          }
          _negativeLine.tcHitsIndices.push_back(j);
          _negativeLine.helixPhiCorrections.push_back(0);
          _negativeLine.zMin = testZ;
          _negativeLine.fitter.addPoint(testZ, testPhi, testWeight);
          lastAddedNegativePhi = testPhi;
          lastAddedNegativePhiError2 = testError2;
        }
      }

      if (_positiveLine.fitter.qn() >= _minLineSegmentHits && _positiveLine.fitter.dydx() <= _maxDphiDz) {
        float dphidz = _positiveLine.fitter.dydx();
        float p2 = computeHelixMomentum2(rC, dphidz);
        if ((seedPhi - lastAddedPositivePhi) * (seedPhi - lastAddedPositivePhi) >
            _segMultiplier * _segMultiplier * seedError2 &&
            p2 > _minHelixMomentum * _minHelixMomentum &&
            p2 < _maxHelixMomentum * _maxHelixMomentum &&
            _positiveLine.fitter.chi2Dof() <= _maxSegmentChi2) {
          _seedPhiLines.push_back(_positiveLine);
          if (_diagLevel == 1) {
            lsInfo.chi2dof = _positiveLine.fitter.chi2Dof();
            lsInfo.maxHitGap = maxHitGap;
            _diagInfo.lineSegmentData.push_back(lsInfo);
          }
          for (size_t n = 0; n < _positiveLine.tcHitsIndices.size(); n++) {
            _tcHits[_positiveLine.tcHitsIndices[n]].notOnSegment = false;
          }
        }
      }
      if (_negativeLine.fitter.qn() >= _minLineSegmentHits && std::abs(_negativeLine.fitter.dydx()) <= _maxDphiDz) {
        float dphidz = _negativeLine.fitter.dydx();
        float p2 = computeHelixMomentum2(rC, dphidz);
        if ((seedPhi - lastAddedNegativePhi) * (seedPhi - lastAddedNegativePhi) >
            _segMultiplier * _segMultiplier * seedError2 &&
            p2 > _minHelixMomentum * _minHelixMomentum &&
            p2 < _maxHelixMomentum * _maxHelixMomentum &&
            _negativeLine.fitter.chi2Dof() <= _maxSegmentChi2) {
          _seedPhiLines.push_back(_negativeLine);
          for (size_t n = 0; n < _negativeLine.tcHitsIndices.size(); n++) {
            _tcHits[_negativeLine.tcHitsIndices[n]].notOnSegment = false;
          }
        }
      }
    }

    // check if we found any lines
    if (_seedPhiLines.size() == 0) { outcome = CONTINUE; }
    else { outcome = GOOD; }

  }

  //-----------------------------------------------------------------------------
  // function to find seed phi lines
  //-----------------------------------------------------------------------------
  float AgnosticHelixFinder::underEstimateSlope(float& phi1, float& phi1Err2, float& phi2, float& phi2Err2, float& dz) {

    float slope = 0.0;
    float phi1Err = std::sqrt(phi1Err2);
    float phi2Err = std::sqrt(phi2Err2);

    if (phi2 > phi1) { slope = std::abs(((phi2-phi2Err)-(phi1+phi1Err))/(dz)); }
    if (phi2 < phi1) { slope = std::abs(((phi2+phi2Err)-(phi1-phi1Err))/(dz)); }

    return slope;

  }

  //-----------------------------------------------------------------------------
  // function to resolve 2 Pi Ambiguities for each seed phi line found
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::resolve2PiAmbiguities() {

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
        int n = std::round(deltaPhi / (2 * M_PI));
        testPhi = testPhi + 2 * M_PI * n;
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
  }

  //-----------------------------------------------------------------------------
  // function for refining the phi line that was found
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::refinePhiLine(size_t lineIndex, bool& removals) {

    float largestResidual2 = 0.0;
    size_t rmIndex = 0;

    for (size_t i = 0; i < _seedPhiLines[lineIndex].tcHitsIndices.size(); i++) {
      size_t tcHitsIndex = _seedPhiLines[lineIndex].tcHitsIndices[i];
      float pointPhi =
        _tcHits[tcHitsIndex].helixPhi + 2 * M_PI * _seedPhiLines[lineIndex].helixPhiCorrections[i];
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

    if (largestResidual2 > _maxPhiZResidual * _maxPhiZResidual) {
      size_t tcHitsIndex = _seedPhiLines[lineIndex].tcHitsIndices[rmIndex];
      float z = getPos(tcHitsIndex).z();
      float phi = _tcHits[tcHitsIndex].helixPhi +
        2 * M_PI * _seedPhiLines[lineIndex].helixPhiCorrections[rmIndex];
      float weight = 1.0 / (_tcHits[tcHitsIndex].helixPhiError2);
      _seedPhiLines[lineIndex].fitter.removePoint(z, phi, weight);
      _seedPhiLines[lineIndex].tcHitsIndices.erase(_seedPhiLines[lineIndex].tcHitsIndices.begin() +
                                                   rmIndex);
      _seedPhiLines[lineIndex].helixPhiCorrections.erase(
                                                         _seedPhiLines[lineIndex].helixPhiCorrections.begin() + rmIndex);
      removals = true;
    }
    else { removals = false; }
  }

  //-----------------------------------------------------------------------------
  // initialize final circle / line seed prior to hit recovery
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::initFinalSeed(LoopCondition& outcome) {

    // make note of number of hits on circle before circle gets updated using the line info
    float seedCircleHits = _circleFitter.qn();

    // first find best phi line
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

    if (_diagLevel == 1) {
      finalLineInfo lInfo;
      lInfo.nHitsRatio = nHitsRatio;
      _diagInfo.lineInfoData.push_back(lInfo);
    }

    // check if we should continue with helix search or not
    if (_circleFitter.qn() < _minFinalSeedHits || nHitsRatio > _maxNHitsRatio) { outcome = CONTINUE;}
    else { outcome = GOOD; }

  }

  //-----------------------------------------------------------------------------
  // recover points
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::recoverPoints(bool& recoveries) {

    float xC = _circleFitter.x0();
    float yC = _circleFitter.y0();
    float rC = _circleFitter.radius();

    float lineSlope = _lineFitter.dydx();
    float lineIntercept = _lineFitter.y0();

    float smallestResidual2 = 0.0;
    size_t addIndex = 0;
    recoveries = false;

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
      int n = std::round(phiDistance / (2 * M_PI));
      _tcHits[i].helixPhiCorrection = n;
      phiP = phiP + 2 * M_PI * n;
      phiDistance = std::abs(lineSlope * zP + lineIntercept - phiP);
      float phiResidual2 = phiDistance * phiDistance / phiSigma2;
      if (phiResidual2 > _maxLineRecoverSigma * _maxLineRecoverSigma) { continue; }
      float residual2 = phiResidual2 + circleResidual2;
      if (residual2 < smallestResidual2 || recoveries == false) {
        smallestResidual2 = residual2;
        addIndex = i;
      }
      recoveries = true;
    }

    // add point to circle fitter to update circle parameters
    // recompute circle points weights
    // update circle with new weights
    // recompute phis / phi errors
    // update line
    if (recoveries == true) {
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
        _tcHits[i].helixPhiCorrection = std::round(deltaPhi / (2 * M_PI));
        phiP = phiP + _tcHits[i].helixPhiCorrection * 2 * M_PI;
        _lineFitter.addPoint(zP, phiP, phiWeight);
      }
    }
  }

  //-----------------------------------------------------------------------------
  // check if helix is savable based on fcl defined cuts
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::checkHelixViability(LoopCondition& outcome) {

    // first check if line has good enough fit
    if (_lineFitter.chi2Dof() > _chi2LineSaveThresh) {
      outcome = CONTINUE;
      return;
    }

    // if line fit was good enough make sure momenta are within allowable range
    float radius = _circleFitter.radius();
    float slope = _lineFitter.dydx();
    float mom2 = computeHelixMomentum2(radius, slope);
    float pt = computeHelixPerpMomentum(radius);
    if (mom2 < _minHelixMomentum * _minHelixMomentum ||
        mom2 > _maxHelixMomentum * _maxHelixMomentum || pt < _minHelixPerpMomentum ||
        pt > _maxHelixPerpMomentum) { outcome = CONTINUE; }
    else { outcome = GOOD; }

  }

  //-----------------------------------------------------------------------------
  // function to save helix
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::saveHelix(size_t tc, HelixSeedCollection& HSColl) {
    HelixSeed hseed;
    hseed._t0 = _tcColl->at(tc)._t0;
    auto _tcCollH = _event->getValidHandle<TimeClusterCollection>(_tcLabel);
    hseed._timeCluster = art::Ptr<mu2e::TimeCluster>(_tcCollH, tc);
    hseed._hhits.setParent(_chColl->parent());

    // flag hits used in helix, and push to combo hit collection in helix seed
    // also add points to linear fitter to get t0
    ::LsqSums2 fitter;
    for (size_t i = 0; i < _tcHits.size(); i++) {
      if (_tcHits[i].inHelix == true || _tcHits[i].hitIndice < 0) { continue; }
      if (_tcHits[i].used == false) { continue; }
      if (_tcHits[i].used == true) { _tcHits[i].inHelix = true; }
      int hitIndice = _tcHits[i].hitIndice;
      const ComboHit* hit = &_chColl->at(hitIndice);
      fitter.addPoint(hit->pos().z(), hit->correctedTime(), 1 / (hit->timeRes() * hit->timeRes()));
      ComboHit hhit(*hit);
      hhit._hphi = _tcHits[i].helixPhi + _tcHits[i].helixPhiCorrection * 2 * M_PI;
      hseed._hhits.push_back(hhit);
    }
    float eDepAvg = hseed._hhits.eDepAvg();

    hseed._t0 = TrkT0(fitter.y0(), fitter.y0Err());

    // set the geom helix parameters
    float x0 = _circleFitter.x0();
    float y0 = _circleFitter.y0();
    CLHEP::Hep3Vector center(x0, y0, 0);
    float dfdz = _lineFitter.dydx();

    hseed._helix._fz0 = _lineFitter.y0();
    if (hseed._helix._fz0 > M_PI) { hseed._helix._fz0 = hseed._helix._fz0 - 2 * M_PI; }
    if (hseed._helix._fz0 < -M_PI) { hseed._helix._fz0 = hseed._helix._fz0 + 2 * M_PI; }
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

    // take care of plotting if _diagLevel = 1
    if (_diagLevel == 1) {
      hsInfo hsi;
      hsi.eDepAvg = eDepAvg;
      _diagInfo.helixSeedData.push_back(hsi);
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
