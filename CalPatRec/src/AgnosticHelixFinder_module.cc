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
#include "Offline/RecoDataProducts/inc/TrkFitFlag.hh"

#include "Offline/Mu2eUtilities/inc/LsqSums2.hh"
#include "Offline/Mu2eUtilities/inc/LsqSums4.hh"
#include "Offline/Mu2eUtilities/inc/polyAtan2.hh"

#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/Mu2eUtilities/inc/McUtilsToolBase.hh"
#include "Offline/Mu2eUtilities/inc/ModuleHistToolBase.hh"

#include "CLHEP/Units/PhysicalConstants.h"

#include "TAxis.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TLine.h"
#include "TMath.h"
#include "TMultiGraph.h"
#include <TROOT.h>

#include <chrono>

namespace mu2e {

  using namespace AgnosticHelixFinderTypes;

  class AgnosticHelixFinder : public art::EDProducer {

  public:
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<int>           diagLevel              {Name("diagLevel"            ), Comment("turn tool on or off"         )  };
      fhicl::Atom<int>           debug                  {Name("debug"                ), Comment("preselect good TC using MC"  )  };
      fhicl::Atom<int>           runDisplay             {Name("runDisplay"           ), Comment("show reco plots when debug"  )  };
      fhicl::Atom<art::InputTag> chCollLabel            {Name("chCollLabel"          ), Comment("combo hit collection label"  )  };
      fhicl::Atom<art::InputTag> tcCollLabel            {Name("tcCollLabel"          ), Comment("time cluster coll label"     )  };
      fhicl::Atom<art::InputTag> ccCollLabel            {Name("ccCollLabel"          ), Comment("Calo Cluster coll label"     )  };
      fhicl::Atom<bool>          findMultipleHelices    {Name("findMultipleHelices"  ), Comment("allow more than one helix"   )  };
      fhicl::Atom<bool>          useStoppingTarget      {Name("useStoppingTarget"    ), Comment("allow in triplet candidates" )  };
      fhicl::Atom<int>           intenseEventThresh     {Name("intenseEventThresh"   ), Comment("# of clusters threshold"     )  };
      fhicl::Atom<int>           intenseClusterThresh   {Name("intenseClusterThresh" ), Comment("# of combo hits threshold"   )  };
      fhicl::Atom<bool>          doIsolationFlag        {Name("doIsolationFlag"      ), Comment("to filter out isolated hits" )  };
      fhicl::Atom<float>         isoRad                 {Name("isoRad"               ), Comment("for isolation cut"           )  };
      fhicl::Atom<int>           isoMinHitsNear         {Name("isoMinHitsNear"       ), Comment("#hits threshold for iso cut" )  };
      fhicl::Atom<bool>          doAverageFlag          {Name("doAverageFlag"        ), Comment("to average out hits or not"  )  };
      fhicl::Atom<float>         minDistCut             {Name("minDistCut"           ), Comment("for averaging out points"    )  };
      fhicl::Atom<float>         minTripletSeedZ        {Name("minTripletSeedZ"      ), Comment("minimum z for triplet seed"  )  };
      fhicl::Atom<float>         minTripletDz           {Name("minTripletDz"         ), Comment("min Z dist btwn 2 trip pts"  )  };
      fhicl::Atom<float>         maxTripletDz           {Name("maxTripletDz"         ), Comment("max Z dist btwn 2 trip pts"  )  };
      fhicl::Atom<float>         minTripletDist         {Name("minTripletDist"       ), Comment("min XY dist btwn 2 trip pts" )  };
      fhicl::Atom<float>         minTripletArea         {Name("minTripletArea"       ), Comment("triangle area of triplet"    )  };
      fhicl::Atom<float>         maxSeedCircleResidual  {Name("maxSeedCircleResidual"), Comment("add hits to triplet circle"  )  };
      fhicl::Atom<int>           minSeedCircleHits      {Name("minSeedCircleHits"    ), Comment("min hits to continue search" )  };
      fhicl::Atom<float>         maxDphiDz              {Name("maxDphiDz"            ), Comment("used finding phi-z segment"  )  };
      fhicl::Atom<float>         maxSeedLineGap         {Name("maxSeedLineGap"       ), Comment("used finding phi-z segment"  )  };
      fhicl::Atom<float>         maxZWindow             {Name("maxZWindow"           ), Comment("used finding phi-z segment"  )  };
      fhicl::Atom<int>           minLineSegmentHits     {Name("minLineSegmentHits"   ), Comment("used in findSeedPhiLines()"  )  };
      fhicl::Atom<float>         segMultiplier          {Name("segMultiplier"        ), Comment("used in findSeedPhiLines()"  )  };
      fhicl::Atom<float>         maxSegmentChi2         {Name("maxSegmentChi2"       ), Comment("used in findSeedPhiLines()"  )  };
      fhicl::Atom<float>         max2PiAmbigResidual    {Name("max2PiAmbigResidual"  ), Comment("when 2pi resolving segment"  )  };
      fhicl::Atom<float>         maxPhiZResidual        {Name("maxPhiZResidual"      ), Comment("when refining phi-z line"    )  };
      fhicl::Atom<int>           minFinalSeedHits       {Name("minFinalSeedHits"     ), Comment("halt search if below thresh" )  };
      fhicl::Atom<float>         maxCircleRecoverSigma  {Name("maxCircleRecoverSigma"), Comment("when doing final recovery"   )  };
      fhicl::Atom<float>         maxLineRecoverSigma    {Name("maxLineRecoverSigma"  ), Comment("when doing final recovery"   )  };
      fhicl::Atom<float>         caloClusterSigma       {Name("caloClusterSigma"     ), Comment("error assigned to calo clust")  };
      fhicl::Atom<int>           minNHelixStrawHits     {Name("minNHelixStrawHits"   ),Comment("straw hit save threshold"     )  };
      fhicl::Atom<int>           minNHelixComboHits     {Name("minNHelixComboHits"   ), Comment("combo hit save threshold"    )  };
      fhicl::Atom<float>         minHelixPerpMomentum   {Name("minHelixPerpMomentum" ), Comment("min pt of helix"             )  };
      fhicl::Atom<float>         maxHelixPerpMomentum   {Name("maxHelixPerpMomentum" ), Comment("max pt of helix"             )  };
      fhicl::Atom<float>         minHelixMomentum       {Name("minHelixMomentum"     ), Comment("min momentum of helix"       )  };
      fhicl::Atom<float>         maxHelixMomentum       {Name("maxHelixMomentum"     ), Comment("max momentum of helix"       )  };
      fhicl::Atom<float>         chi2LineSaveThresh     {Name("chi2LineSaveThresh"   ), Comment("max chi2Dof for line"        )  };
      fhicl::Atom<float>         maxEDepAvg             {Name("maxEDepAvg"           ), Comment("max avg edep of combohits"   )  };
      fhicl::Atom<float>         maxNHitsRatio          {Name("maxNHitsRatio"        ), Comment("max ratio of seed hits"      )  };
      fhicl::Atom<int>           debugPdgID             {Name("debugPdgID"           ), Comment("pdgID of interest in display")  };
      fhicl::Atom<float>         debugMomentum          {Name("debugMomentum"        ), Comment("lower momentum of interest"  )  };
      fhicl::Atom<std::string>   debugDirection         {Name("debugDirection"       ), Comment("down or up"                  )  };
      fhicl::Atom<int>           debugStrawHitThresh    {Name("debugStrawHitThresh"  ), Comment("min #SHs for plot particle"  )  };
      fhicl::Atom<float>         debugScatterThresh     {Name("debugScatterThresh"   ), Comment("max dP for plot particle"    )  };

      fhicl::Table<McUtilsToolBase::Config>          mcUtils     {Name("mcUtils"   ), Comment("get MC info if debugging"      )  };
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
      bool    notOnLine = true;
      bool    uselessTripletSeed = false;
      bool    notOnSegment = true;
      bool    debugParticle = false; // only filled in debug mode -- true if mc particle, false if background
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

    struct helixCandidate {
      std::vector<cHit>   tcHitsCopy;
      ::LsqSums4          circleFitterCopy;
      ::LsqSums2          lineFitterCopy;
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

    // struct to hold info specifically for when debugging
    struct logicTrackingInfo {
      std::vector<std::vector<cHit>>   tcHitsColl;
      std::vector<::LsqSums4>          circleFitter;
      std::vector<::LsqSums2>          lineFitter;
      std::vector<lineInfo>            lineSegments;
      std::vector<lineInfo>            resolvedSegments;
      bool                             foundLineSegment;
      size_t                           bestLineSegment;
      int                              finalNHits;
      bool                             becameHelix;
    };

    // another struct for debugging
    struct mcInfo {
      int   simID;
      int   nStrawHits;
      float pMin;
      float pMax;
    };

  private:
    //-----------------------------------------------------------------------------
    // tool on or off, debug and runDisplay fcl parameters
    //-----------------------------------------------------------------------------
    int _diagLevel;
    int _debug;
    int _runDisplay;

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
    float    _maxNHitsRatio;

    //-----------------------------------------------------------------------------
    // fcl params to choose what particle to plot when doing debugging
    //-----------------------------------------------------------------------------
    int           _debugPdgID;
    float         _debugMomentum;
    std::string   _debugDirection;
    int           _debugStrawHitThresh;
    float         _debugScatterThresh;

    //-----------------------------------------------------------------------------
    // need to use mcUtils if in debug mode
    //-----------------------------------------------------------------------------
    std::unique_ptr<McUtilsToolBase> _mcUtils;

    //-----------------------------------------------------------------------------
    // diagnostics
    //-----------------------------------------------------------------------------
    art::Handle<CaloClusterCollection>   _ccHandle;
    const mu2e::Calorimeter*             _calorimeter;
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
    std::vector<helixCandidate>   _helixCandidates;
    float                         _bz0;
    bool                          _intenseEvent;
    bool                          _intenseCluster;

    //-----------------------------------------------------------------------------
    // data members specifically for when doing debugging
    //-----------------------------------------------------------------------------
    std::vector<logicTrackingInfo>     _plottingData;
    std::vector<std::vector<mcInfo>>   _simIDsPerTC; // filled once per TC
    int                                _tcIndex;
    int                                _simID;
    float                              _mcRadius;
    float                              _mcX0;
    float                              _mcY0;
    size_t                             _bestPlotIndex;
    size_t                             _bestLineSegment;
    TCanvas*                           _xy0;
    TCanvas*                           _xy1;
    TCanvas*                           _xy2;
    TCanvas*                           _phiz2;
    TCanvas*                           _bestLineSeg;
    TCanvas*                           _bestLineSegResolved;
    TCanvas*                           _xy3;
    TCanvas*                           _phiz3;
    TCanvas*                           _xy4;
    TCanvas*                           _phiz4;

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
    void         findSeedPhiLines          ();
    float        underEstimateSlope        (float& phi1, float& phi1Err2, float& phi2, float& phi2Err2, float& dz);
    void         resolve2PiAmbiguities     ();
    void         refinePhiLine             (size_t lineIndex, bool& removals);
    void         initFinalSeed             ();
    void         recoverPoints             (bool& recoveries);
    void         saveHelixCandidate        ();
    void         saveHelix                 (size_t tc, HelixSeedCollection& HSColl);

    //-----------------------------------------------------------------------------
    // functions for debug mode and runDisplay mode
    //-----------------------------------------------------------------------------
    void initDebugMode      ();
    void findBestTC         ();
    void findBestPlotIndex  ();
    void doAllPlots         ();
    void plotXY             (int stage);
    void plotPhiZ           (int stage);
    void plotSegment        (int option);
  };

  //-----------------------------------------------------------------------------
  // module constructor
  //-----------------------------------------------------------------------------
  AgnosticHelixFinder::AgnosticHelixFinder(const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config},
    _diagLevel                     (config().diagLevel()                             ),
    _debug                         (config().debug()                                 ),
    _runDisplay                    (config().runDisplay()                            ),
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
    _maxNHitsRatio                 (config().maxNHitsRatio()                         ),
    _debugPdgID                    (config().debugPdgID()                            ),
    _debugMomentum                 (config().debugMomentum()                         ),
    _debugDirection                (config().debugDirection()                        ),
    _debugStrawHitThresh           (config().debugStrawHitThresh()                   ),
    _debugScatterThresh            (config().debugScatterThresh()                    )

    {

      consumes<ComboHitCollection>     (_chLabel);
      consumes<TimeClusterCollection>  (_tcLabel);
      consumes<CaloClusterCollection>  (_ccLabel);
      produces<HelixSeedCollection>    ();

      if (_debug == 1) {
        _mcUtils = art::make_tool<McUtilsToolBase>(config().mcUtils, "mcUtils");
        if (_runDisplay == 1) {
          _xy0                 = new TCanvas("_xy0", "XY view of all hits in time cluster (stage 0)", 900, 900);
          _xy1                 = new TCanvas("_xy1", "XY view of triplet (stage 1)", 900, 900);
          _xy2                 = new TCanvas("_xy2", "XY view of seed circle (stage 2)", 900, 900);
          _phiz2               = new TCanvas("_phiz2", "Phi-Z view after computing phi relative to circle center (stage 2)", 900, 900);
          _bestLineSeg         = new TCanvas("_bestLineSeg", "Line segment that leads to best candidate line", 900, 900);
          _bestLineSegResolved = new TCanvas("_bestLineSegResolved","2pi resolved segment that leads to best candidate", 900, 900);
          _xy3                 = new TCanvas("_xy3", "XY view after updating based on line found (stage 3)", 900, 900);
          _phiz3               = new TCanvas("_phiz3", "Phi-Z view of best line candidate found (stage 3)", 900, 900);
          _xy4                 = new TCanvas("_xy4", "XY view of final circle after recovery stage (stage 4)", 900, 900);
          _phiz4               = new TCanvas("_phiz4", "Phi-z view of final line after recovery stage (stage 4)", 900, 900);
        }
      }
      else _mcUtils = std::make_unique<McUtilsToolBase>();

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
    } else {
      _chColl = 0;
    }

    auto _tcCollH = evt.getValidHandle<TimeClusterCollection>(_tcLabel);
    if (_tcCollH.product() != 0) {
      _tcColl = _tcCollH.product();
    } else {
      _tcColl = 0;
    }

    if (evt.getByLabel(_ccLabel, _ccHandle)) {
      _ccColl = _ccHandle.product();
    } else {
      _ccColl = NULL;
    }

    return (_tcColl != 0);
  }

  //-----------------------------------------------------------------------------
  // event entry point
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::produce(art::Event& event) {

    // get time (needed for diagnostic tool)
    auto moduleStartTime = std::chrono::high_resolution_clock::now();

    // set the event and prepare hsColl
    _event = &event;
    std::unique_ptr<HelixSeedCollection> hsColl(new HelixSeedCollection);

    // get the necessary data to do helix search
    bool dataExists = findData(event);

    // prepare diagnostic tool data members
    if (_diagLevel == 1) {
      _diagInfo.moduleTime = 0.0;
      _diagInfo.nHelices = 0;
      _diagInfo.timeClusterData.clear();
      _diagInfo.helixSeedData.clear();
      _diagInfo.lineSegmentData.clear();
    }

    // flag whether event is intense or not
    if ((int)_tcColl->size() > _intenseEventThresh) {
      _intenseEvent = true;
    } else {
      _intenseEvent = false;
    }

    // run necessary logic if in debug mode
    if (_debug == 1) {
      initDebugMode();
    } else {
      _tcIndex = -1;
    }

    // do helix search
    bool continueSearch = true;

    // search for helix
    if (dataExists) {
      for (size_t i = 0; i < _tcColl->size(); i++) {
        if (_debug == 1 && (int)i != _tcIndex) {
          continue;
        } // only search for helix in TC of interest if in debug mode
        // check to see if cluster is a busy one
        _intenseCluster = false;
        if ((int)_tcColl->at(i).nhits() > _intenseClusterThresh) {
          _intenseCluster = true;
        }
        // if cluster is busy or event is intense then only process time cluster if it has calo
        // cluster
        if (_intenseEvent == true || _intenseCluster == true) {
          if (_tcColl->at(i).hasCaloCluster() == false) {
            continue;
          }
        }
        tcInfo timeClusterInfo;
        auto tcStartTime = std::chrono::high_resolution_clock::now();
        int nHelicesInitial = _diagInfo.nHelices;
        _tcHits.clear();
        _helixCandidates.clear();
        tcHitsFill(i);
        continueSearch = true;
        while (continueSearch == true) {
          findHelix(i, *hsColl, continueSearch);
          if (_findMultipleHelices == false) {
            continueSearch = false;
          } // want to halt search if not configured to find multiple helices per TC
        }
        if (_diagLevel == 1) {
          auto tcEndTime = std::chrono::high_resolution_clock::now();
          timeClusterInfo.time =
            std::chrono::duration<float, std::milli>(tcEndTime - tcStartTime).count();
          timeClusterInfo.nHelices = _diagInfo.nHelices - nHelicesInitial;
          timeClusterInfo.nComboHits = _tcColl->at(i).nhits();
          timeClusterInfo.nStrawHits = _tcColl->at(i).nStrawHits();
          _diagInfo.timeClusterData.push_back(timeClusterInfo);
        }
      }
      if (_runDisplay == 1 && _tcIndex != -1) {
        doAllPlots();
      } // if in runDisplay mode and we have good TC, then make display plots
    }

    // put helix seed collection into the event record
    event.put(std::move(hsColl));

    // fill necessary data members for diagnostic tool
    if (_diagLevel == 1) {
      auto moduleEndTime = std::chrono::high_resolution_clock::now();
      _diagInfo.moduleTime =
        std::chrono::duration<float, std::milli>(moduleEndTime - moduleStartTime).count();
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

    if (hitIndice >= 0) {
      return _chColl->at(hitIndice).pos();
    }
    if (hitIndice == HitType::STOPPINGTARGET) {
      return _stopTargPos;
    }
    if (hitIndice == HitType::CALOCLUSTER) {
      return _caloPos;
    }

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
        _chColl->at(hitIndice).wireVar() * sinth2 + transVar * costh2;
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
      if (_debug == 1) {
        std::vector<StrawDigiIndex> shids;
        _chColl->fillStrawDigiIndices(hit.hitIndice, shids);
        int SimID = _mcUtils->strawHitSimId(_event, shids[0]);
        if (SimID == _simID) {
          hit.debugParticle = true;
        } else {
          hit.debugParticle = false;
        }
      }
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

    // do isolation and average flagging
    if (_doIsolationFlag == true || _doAverageFlag == true) {
      for (size_t i = 0; i < _tcHits.size(); i++) {
        if (_tcHits[i].inHelix == true || _tcHits[i].hitIndice < 0) {
          continue;
        }
        int nHitsNear = 0;
        XYZVectorF seedPos = getPos(i);
        for (size_t j = 0; j < _tcHits.size(); j++) {
          if (_tcHits[j].inHelix == true || _tcHits[j].hitIndice < 0) {
            continue;
          }
          if (j == i) {
            continue;
          }
          XYZVectorF testPos = getPos(j);
          // do isolation flagging
          if (_doIsolationFlag == true) {
            if ((seedPos-testPos).Perp2() < _isoRad * _isoRad) {
              nHitsNear++;
            }
            if (nHitsNear < _isoMinHitsNear) {
              _tcHits[i].isolated = true;
            } else {
              _tcHits[i].isolated = false;
            }
          }
          // do averaging out
          if (_doAverageFlag == true) {
            if (_tcHits[i].averagedOut == true) {
              continue;
            }
            if (_tcHits[j].averagedOut == true) {
              continue;
            }
            if ((seedPos-testPos).Perp2() <= _minDistCut * _minDistCut) {
              _tcHits[j].averagedOut = true;
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

    // first flag hits that are isolated
    for (size_t i = 0; i < _tcHits.size(); i++) {
      if (_tcHits[i].inHelix == true || _tcHits[i].hitIndice == HitType::STOPPINGTARGET) {
        continue;
      }
      _tcHits[i].isolated = false;
      _tcHits[i].averagedOut = false;
      _tcHits[i].notOnLine = true;
      _tcHits[i].notOnSegment = true;
    }
  }

  //-----------------------------------------------------------------------------
  // logic to find helix
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::findHelix(size_t tc, HelixSeedCollection& HSColl, bool& findAnotherHelix) {

    // if in runDisplay mode we store collections for plotting
    _circleFitter.clear();
    _lineFitter.clear();
    logicTrackingInfo trackingInfo;
    if (_debug == 1 && _runDisplay == 1) { // stage 0 info (all hits in TC)
      trackingInfo.tcHitsColl.push_back(_tcHits);
      trackingInfo.circleFitter.push_back(_circleFitter);
      trackingInfo.lineFitter.push_back(_lineFitter);
    }

    // by default we set findAnotherHelix to false, if we find a helix then this will be set to true
    findAnotherHelix = false;

    // set flags before starting search so that we know what hits to use for tripletting
    if (_doIsolationFlag == true || _doAverageFlag == true) {
      setFlags();
    }

    // now we loop over triplets
    bool foundPhiZRemoval = true;
    bool pointRecovered = true;
    for (size_t i = 0; i < _tcHits.size() - 2; i++) {
      bool uselessSeed = true;
      triplet tripletInfo;
      LoopCondition loopCondition;
      setTripletI(i, tripletInfo, loopCondition);
      if (loopCondition == CONTINUE) {
        continue;
      }
      if (loopCondition == BREAK) {
        break;
      }
      for (size_t j = i + 1; j < _tcHits.size() - 1; j++) {
        setTripletJ(j, tripletInfo, loopCondition);
        if (loopCondition == CONTINUE) {
          continue;
        }
        if (loopCondition == BREAK) {
          break;
        }
        for (size_t k = j + 1; k < _tcHits.size(); k++) {
          setTripletK(k, tripletInfo, loopCondition);
          if (loopCondition == CONTINUE) {
            continue;
          }
          if (loopCondition == BREAK) {
            break;
          }
          // clear fitters for new triplet then initialize triplet
          _circleFitter.clear();
          _lineFitter.clear();
          initTriplet(tripletInfo, loopCondition);
          // if in debug / runDisplay mode we want to make plot with triplet
          if (_debug == 1 && _runDisplay == 1) { // stage 1 info (triplet)
            for (size_t n = 0; n < _tcHits.size(); n++) {
              if (n == i) {
                _tcHits[n].used = true;
              }
              if (n == j) {
                _tcHits[n].used = true;
              }
              if (n == k) {
                _tcHits[n].used = true;
              }
              if (n != i && n != j && n != k) {
                _tcHits[n].used = false;
              }
            }
            trackingInfo.tcHitsColl.push_back(_tcHits);
            trackingInfo.circleFitter.push_back(_circleFitter);
            trackingInfo.lineFitter.push_back(_lineFitter);
          }
          // now initialize seed circle if triplet circle passed condition check
          if (loopCondition == CONTINUE) {
            continue;
          }
          initSeedCircle(loopCondition);
          float circleHits = _circleFitter.qn();
          if (loopCondition == CONTINUE) {
            if (_debug == 1 && _runDisplay == 1) { // stage 2 info (seed circle + helix phi)
              initHelixPhi();
              trackingInfo.tcHitsColl.push_back(_tcHits);
              trackingInfo.circleFitter.push_back(_circleFitter);
              trackingInfo.lineFitter.push_back(_lineFitter);
              trackingInfo.finalNHits = _circleFitter.qn();
              trackingInfo.becameHelix = false;
              _plottingData.push_back(trackingInfo);
            }
            continue;
          }
          uselessSeed = false;
          initHelixPhi();
          if (_debug == 1 && _runDisplay == 1) { // stage 2 info (seed circle + helix phi)
            trackingInfo.tcHitsColl.push_back(_tcHits);
            trackingInfo.circleFitter.push_back(_circleFitter);
            trackingInfo.lineFitter.push_back(_lineFitter);
          }
          findSeedPhiLines();
          if (_seedPhiLines.size() == 0) {
            if (_debug == 1 && _runDisplay == 1) {
              trackingInfo.foundLineSegment = false;
              trackingInfo.finalNHits = _circleFitter.qn();
              trackingInfo.becameHelix = false;
              _plottingData.push_back(trackingInfo);
            }
            continue;
          }
          if (_debug == 1 && _runDisplay == 1) {
            trackingInfo.foundLineSegment = true;
            trackingInfo.lineSegments = _seedPhiLines;
          }
          resolve2PiAmbiguities();
          if (_debug == 1 && _runDisplay == 1) {
            trackingInfo.resolvedSegments = _seedPhiLines;
          }
          for (size_t i = 0; i < _seedPhiLines.size(); i++) {
            if ((int)_seedPhiLines[i].tcHitsIndices.size() < _minFinalSeedHits) {
              continue;
            }
            foundPhiZRemoval = true;
            while (foundPhiZRemoval == true) {
              refinePhiLine(i, foundPhiZRemoval);
            }

          }

          initFinalSeed();
          float phiHits = _lineFitter.qn();
          float nHitsRatio = circleHits/(phiHits + 1e-6);
          if (_debug == 1) {
            hsInfo helixSeedInfo;
            helixSeedInfo.nHitsRatio = nHitsRatio;
            _diagInfo.helixSeedData.push_back(helixSeedInfo);
          }
          if (_debug == 1 && _runDisplay == 1) { // stage 3 info (final seed)
            trackingInfo.bestLineSegment = _bestLineSegment;
            trackingInfo.tcHitsColl.push_back(_tcHits);
            trackingInfo.circleFitter.push_back(_circleFitter);
            trackingInfo.lineFitter.push_back(_lineFitter);
          }
          if (_circleFitter.qn() < _minFinalSeedHits) {
            if (_debug == 1 && _runDisplay == 1) {
              trackingInfo.finalNHits = _circleFitter.qn();
              trackingInfo.becameHelix = false;
              _plottingData.push_back(trackingInfo);
            }
            continue;
          }
          if (nHitsRatio > _maxNHitsRatio) {
            continue;
          }
          pointRecovered = true;
          while (pointRecovered == true) {
            recoverPoints(pointRecovered);
          }
          if (_debug == 1 && _runDisplay == 1) { // stage 4 info (final helix)
            trackingInfo.tcHitsColl.push_back(_tcHits);
            trackingInfo.circleFitter.push_back(_circleFitter);
            trackingInfo.lineFitter.push_back(_lineFitter);
          }
          if (_lineFitter.chi2Dof() > _chi2LineSaveThresh) {
            if (_debug == 1 && _runDisplay == 1) {
              trackingInfo.finalNHits = _circleFitter.qn();
              trackingInfo.becameHelix = false;
              _plottingData.push_back(trackingInfo);
            }
            continue;
          }
          // see if helix is within momentum range allowed
          float radius = _circleFitter.radius();
          float slope = _lineFitter.dydx();
          float mom2 = computeHelixMomentum2(radius, slope);
          float pt = computeHelixPerpMomentum(radius);
          if (mom2 < _minHelixMomentum * _minHelixMomentum ||
              mom2 > _maxHelixMomentum * _maxHelixMomentum || pt < _minHelixPerpMomentum ||
              pt > _maxHelixPerpMomentum) {
            if (_debug == 1 && _runDisplay == 1) {
              trackingInfo.finalNHits = _circleFitter.qn();
              trackingInfo.becameHelix = false;
              _plottingData.push_back(trackingInfo);
            }
            continue;
          }
          int nStrawHitsInHelix = 0;
          int nComboHitsInHelix = 0;
          int nStrawHitsInTimeCluster = 0;
          int nComboHitsInTimeCluster = 0;
          for (size_t q = 0; q < _tcHits.size(); q++) {
            if (_tcHits[q].inHelix == true || _tcHits[q].hitIndice < 0) {
              continue;
            }
            int hitIndice = _tcHits[q].hitIndice;
            nStrawHitsInTimeCluster = nStrawHitsInTimeCluster + _chColl->at(hitIndice).nStrawHits();
            nComboHitsInTimeCluster = nComboHitsInTimeCluster + 1;
            if (_tcHits[q].used == false) {
              continue;
            }
            nStrawHitsInHelix = nStrawHitsInHelix + _chColl->at(hitIndice).nStrawHits();
            nComboHitsInHelix = nComboHitsInHelix + 1;
          }
          if (nStrawHitsInHelix >= _minNHelixStrawHits && nComboHitsInHelix >= _minNHelixComboHits) {
            if (_debug == 1 && _runDisplay == 1) {
              trackingInfo.finalNHits = _circleFitter.qn();
              trackingInfo.becameHelix = true;
              _plottingData.push_back(trackingInfo);
            }
            saveHelixCandidate();
            saveHelix(tc, HSColl);
            if (_diagLevel == 1) {
              _diagInfo.nHelices++;
            }
            int remainingStrawHits = nStrawHitsInTimeCluster - nStrawHitsInHelix;
            int remainingComboHits = nComboHitsInTimeCluster - nComboHitsInHelix;
            if (remainingStrawHits < _minNHelixStrawHits ||
                remainingComboHits < _minNHelixComboHits) {
              i = _tcHits.size();
              j = i;
              break;
            } else {
              findAnotherHelix = true;
              if (_doIsolationFlag == true || _doAverageFlag == true) {
                resetFlags();
              }
              i = _tcHits.size();
              j = i;
              break;
            }
          } else {
            if (_debug == 1 && _runDisplay == 1) {
              trackingInfo.finalNHits = _circleFitter.qn();
              trackingInfo.becameHelix = false;
              _plottingData.push_back(trackingInfo);
            }
          }
        }
      }
      if (i<_tcHits.size()) {_tcHits[i].uselessTripletSeed = uselessSeed;}
    }
  }

  //-----------------------------------------------------------------------------
  // check flags to see if point is good for triplet-ing with
  //-----------------------------------------------------------------------------
  bool AgnosticHelixFinder::passesFlags(size_t& tcHitsIndex) {

    if (_tcHits[tcHitsIndex].inHelix == true) {
      return false;
    }
    if (_tcHits[tcHitsIndex].uselessTripletSeed == true) {
      return false;
    }
    if (_doIsolationFlag == true && _tcHits[tcHitsIndex].isolated == true) {
      return false;
    }
    if (_doAverageFlag == true && _tcHits[tcHitsIndex].averagedOut == true) {
      return false;
    }

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
    if (pt < _minHelixPerpMomentum || pt > _maxHelixPerpMomentum) {
      outcome = CONTINUE;
    } else {
      outcome = GOOD;
    }
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
      if (_tcHits[i].inHelix == true) {
        continue;
      }
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
    if (_circleFitter.qn() < _minSeedCircleHits) {
      outcome = CONTINUE;
    } else {
      outcome = GOOD;
    }
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
  void AgnosticHelixFinder::findSeedPhiLines() {

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
      if (_tcHits[i].inHelix == true || _tcHits[i].hitIndice < 0) {
        continue;
      }
      if (_tcHits[i].used == false) {
        continue;
      }
      if (_tcHits[i].notOnSegment == false) {
        continue;
      }
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
        if (_tcHits[j].inHelix == true || _tcHits[j].hitIndice < 0) {
          continue;
        }
        if (_tcHits[j].used == false) {
          continue;
        }
        float testZ = getPos(j).z();
        if (testZ == seedZ) {
          continue;
        }
        if (std::abs(seedZ - testZ) > _maxZWindow) {
          break;
        }
        float positiveDeltaZ = lastAddedPositiveZ - testZ;
        float negativeDeltaZ = lastAddedNegativeZ - testZ;
        if (positiveDeltaZ > _maxSeedLineGap && negativeDeltaZ > _maxSeedLineGap) {
          break;
        }
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
          if (positiveDeltaZ > maxHitGap) {
            maxHitGap = positiveDeltaZ;
          }
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
        if (alreadyInLine == true) {
          continue;
        }
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
          if (testZ < _seedPhiLines[i].zMin) {
            _seedPhiLines[i].zMin = testZ;
          }
          if (testZ > _seedPhiLines[i].zMax) {
            _seedPhiLines[i].zMax = testZ;
          }
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
    } else {
      removals = false;
    }
  }

  //-----------------------------------------------------------------------------
  // initialize final circle / line seed prior to hit recovery
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::initFinalSeed() {

    // first find best phi line
    size_t bestLineIndex = 0;
    size_t nHitsMost = 0;
    float chi2dof = 0.0;
    for (size_t i = 0; i < _seedPhiLines.size(); i++) {
      if (_seedPhiLines[i].tcHitsIndices.size() < nHitsMost) {
        continue;
      }
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

    if (_debug == 1 && _runDisplay == 1) {
      _bestLineSegment = bestLineIndex;
    }

    // now initialize line/circle fitters and flag hits in tcHits that are used
    _circleFitter.clear();
    _lineFitter.clear();
    _lineFitter = _seedPhiLines[bestLineIndex].fitter;
    for (size_t i = 0; i < _tcHits.size(); i++) {
      if (_tcHits[i].inHelix == true || _tcHits[i].hitIndice == HitType::STOPPINGTARGET) {
        continue;
      }
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
      if (phiResidual2 > _maxLineRecoverSigma * _maxLineRecoverSigma) {
        continue;
      }
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
  // function to make copies of relevant info for helix candidate
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::saveHelixCandidate() {

    helixCandidate hel;
    hel.tcHitsCopy = _tcHits;
    hel.lineFitterCopy = _lineFitter;
    hel.circleFitterCopy = _circleFitter;

    _helixCandidates.push_back(hel);
  }

  //-----------------------------------------------------------------------------
  // function to save helix
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::saveHelix(size_t tc, HelixSeedCollection& HSColl) {

    size_t bestIndex = 0;
    _tcHits = _helixCandidates[bestIndex].tcHitsCopy;
    _lineFitter = _helixCandidates[bestIndex].lineFitterCopy;
    _circleFitter = _helixCandidates[bestIndex].circleFitterCopy;
    _helixCandidates.clear();

    HelixSeed hseed;
    hseed._t0 = _tcColl->at(tc)._t0;
    auto _tcCollH = _event->getValidHandle<TimeClusterCollection>(_tcLabel);
    hseed._timeCluster = art::Ptr<mu2e::TimeCluster>(_tcCollH, tc);
    hseed._hhits.setParent(_chColl->parent());

    float eDepSum = 0.0;
    size_t nComboHits = 0;

    // flag hits used in helix, and push to combo hit collection in helix seed
    // also add points to linear fitter to get t0
    ::LsqSums2 fitter;
    for (size_t i = 0; i < _tcHits.size(); i++) {
      if (_tcHits[i].inHelix == true || _tcHits[i].hitIndice < 0) {
        continue;
      }
      if (_tcHits[i].used == false) {
        continue;
      }
      if (_tcHits[i].used == true) {
        _tcHits[i].inHelix = true;
      }
      int hitIndice = _tcHits[i].hitIndice;
      const ComboHit* hit = &_chColl->at(hitIndice);
      fitter.addPoint(hit->pos().z(), hit->correctedTime(), 1 / (hit->timeRes() * hit->timeRes()));
      ComboHit hhit(*hit);
      hhit._hphi = _tcHits[i].helixPhi + _tcHits[i].helixPhiCorrection * 2 * M_PI;
      hseed._hhits.push_back(hhit);
      eDepSum = eDepSum + hhit.energyDep();
      nComboHits++;
    }
    float eDepAvg = eDepSum/nComboHits;

    hseed._t0 = TrkT0(fitter.y0(), fitter.y0Err());

    // set the geom helix parameters
    float x0 = _circleFitter.x0();
    float y0 = _circleFitter.y0();
    CLHEP::Hep3Vector center(x0, y0, 0);
    float dfdz = _lineFitter.dydx();

    hseed._helix._fz0 = _lineFitter.y0();
    if (hseed._helix._fz0 > M_PI) {
      hseed._helix._fz0 = hseed._helix._fz0 - 2 * M_PI;
    }
    if (hseed._helix._fz0 < -M_PI) {
      hseed._helix._fz0 = hseed._helix._fz0 + 2 * M_PI;
    }
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
      _diagInfo.helixSeedData[bestIndex].eDepAvg = eDepAvg;
    }

    if (eDepAvg > _maxEDepAvg) return;

    // push back the helix seed to the helix seed collection
    HSColl.emplace_back(hseed);

  }

  //-----------------------------------------------------------------------------
  // calling logic that needs to be called to run debug mode
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::initDebugMode() {

    // first find the TC we want to focus on
    findBestTC();

    // now prepare display plots if in display mode
    if (_runDisplay == 1) {
      _plottingData.clear();
      _xy0->Destructor();
      _xy0 = new TCanvas("_xy0", "XY view of all hits in time cluster (stage 0)", 900, 900);
      _xy0->SetLeftMargin(0.145);
      _xy1->Destructor();
      _xy1 = new TCanvas("_xy1", "XY view of triplet (stage 1)", 900, 900);
      _xy1->SetLeftMargin(0.145);
      _xy2->Destructor();
      _xy2 = new TCanvas("_xy2", "XY view of seed circle (stage 2)", 900, 900);
      _xy2->SetLeftMargin(0.145);
      _phiz2->Destructor();
      _phiz2 = new TCanvas(
                           "_phiz5", "Phi-Z view after computing phi relative to circle center (stage 2)", 900, 900);
      _phiz2->SetLeftMargin(0.145);
      _bestLineSeg->Destructor();
      _bestLineSeg = new TCanvas(
                                 "_bestLineSeg", "Phi-Z view of line segment that leads to best candidate line", 900, 900);
      _bestLineSeg->SetLeftMargin(0.145);
      _bestLineSegResolved->Destructor();
      _bestLineSegResolved = new TCanvas(
                                         "_bestLineSegResolved",
                                         "Phi-Z view of 2pi resolved line segment that leads to best candidate line", 900, 900);
      _bestLineSegResolved->SetLeftMargin(0.145);
      _xy3->Destructor();
      _xy3 = new TCanvas("_xy3", "XY view after updating based on line found (stage 3)", 900, 900);
      _xy3->SetLeftMargin(0.145);
      _phiz3->Destructor();
      _phiz3 = new TCanvas("_phiz3", "Phi-Z view of best line candidate found (stage 3)", 900, 900);
      _phiz3->SetLeftMargin(0.145);
      _xy4->Destructor();
      _xy4 = new TCanvas("_xy4", "XY view of final circle after recovery stage (stage 4)", 900, 900);
      _xy4->SetLeftMargin(0.145);
      _phiz4->Destructor();
      _phiz4 =
        new TCanvas("_phiz4", "Phi-z view of final line after recovery stage (stage 4)", 900, 900);
      _phiz4->SetLeftMargin(0.145);
    }
  }

  //-----------------------------------------------------------------------------
  // function to find the best time cluster in debug mode given particle of interest (set in fcl)
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::findBestTC() {

    _simIDsPerTC.clear();

    GlobalConstantsHandle<ParticleDataList> pdt;
    float qSign = pdt->particle(_debugPdgID).charge();

    const mu2e::SimParticle* _simParticle;
    // loop over TCs to fill _simIDsPerTC
    for (size_t i = 0; i < _tcColl->size(); i++) {
      std::vector<mcInfo> particlesInTC;
      for (size_t j = 0; j < _tcColl->at(i)._strawHitIdxs.size(); j++) {
        int hitIndice = _tcColl->at(i)._strawHitIdxs[j];
        std::vector<StrawDigiIndex> shids;
        _chColl->fillStrawDigiIndices(hitIndice, shids);
        for (size_t k = 0; k < shids.size(); k++) {
          _simParticle = _mcUtils->getSimParticle(_event, shids[k]);
          if (_simParticle->pdgId() != _debugPdgID) {
            continue;
          }
          const XYZVectorF* simMomentum = _mcUtils->getMom(_event, shids[k]);
          float simXmomentum = simMomentum->x();
          float simYmomentum = simMomentum->y();
          float simZmomentum = simMomentum->z();
          if (_debugDirection == "down" && simZmomentum < 0.0) {
            continue;
          }
          if (_debugDirection == "up" && simZmomentum > 0.0) {
            continue;
          }
          float simPerpMomentum =
            std::sqrt(simXmomentum * simXmomentum + simYmomentum * simYmomentum);
          float simMomentumMag =
            std::sqrt(simPerpMomentum * simPerpMomentum + simZmomentum * simZmomentum);
          if (simMomentumMag < _debugMomentum) {
            continue;
          }
          int SimID = _mcUtils->strawHitSimId(_event, shids[k]);
          bool particleAlreadyFound = false;
          for (size_t n = 0; n < particlesInTC.size(); n++) {
            if (SimID == particlesInTC[n].simID) {
              particleAlreadyFound = true;
              particlesInTC[n].nStrawHits = particlesInTC[n].nStrawHits + 1;
              if (simMomentumMag > particlesInTC[n].pMax) {
                particlesInTC[n].pMax = simMomentumMag;
              }
              if (simMomentumMag < particlesInTC[n].pMin) {
                particlesInTC[n].pMin = simMomentumMag;
              }
              break;
            }
          }
          if (particleAlreadyFound) {
            continue;
          }
          mcInfo particle;
          particle.simID = SimID;
          particle.nStrawHits = 1;
          particle.pMax = simMomentumMag;
          particle.pMin = simMomentumMag;
          particlesInTC.push_back(particle);
        }
      }
      _simIDsPerTC.push_back(particlesInTC);
    }

    // loop over _simIDsPerTC to find best TC
    _tcIndex = -1;
    int mostStrawHits = 0;
    for (size_t i = 0; i < _simIDsPerTC.size(); i++) {
      for (size_t j = 0; j < _simIDsPerTC[i].size(); j++) {
        if (_simIDsPerTC[i].at(j).nStrawHits < _debugStrawHitThresh) {
          continue;
        }
        float momentumDiff = _simIDsPerTC[i].at(j).pMax - _simIDsPerTC[i].at(j).pMin;
        if (momentumDiff > _debugScatterThresh) {
          continue;
        }
        if (_simIDsPerTC[i].at(j).nStrawHits > mostStrawHits) {
          mostStrawHits = _simIDsPerTC[i].at(j).nStrawHits;
          _simID = _simIDsPerTC[i].at(j).simID;
          _tcIndex = (int)i;
        }
      }
    }

    // compute helix MC circle parameters
    if (_tcIndex != -1) {
      for (size_t j = 0; j < _tcColl->at(_tcIndex)._strawHitIdxs.size(); j++) {
        int hitIndice = _tcColl->at(_tcIndex)._strawHitIdxs[j];
        std::vector<StrawDigiIndex> shids;
        _chColl->fillStrawDigiIndices(hitIndice, shids);
        bool foundParticle = false;
        for (size_t k = 0; k < shids.size(); k++) {
          if (_mcUtils->strawHitSimId(_event, shids[k]) != _simID) {
            continue;
          } else {
            foundParticle = true;
            const XYZVectorF* simMomentum = _mcUtils->getMom(_event, shids[k]);
            const XYZVectorF* simPosition = _mcUtils->getPos(_event, shids[k]);
            float simXmomentum = simMomentum->x();
            float simYmomentum = simMomentum->y();
            float simPerpMomentum =
              std::sqrt(simXmomentum * simXmomentum + simYmomentum * simYmomentum);
            _mcRadius = (simPerpMomentum) / (_bz0 * mmTconversion);
            float hitX = simPosition->x();
            float hitY = simPosition->y();
            ;
            _mcX0 = hitX + qSign * simYmomentum * _mcRadius / simPerpMomentum;
            _mcY0 = hitY - qSign * simXmomentum * _mcRadius / simPerpMomentum;
            break;
          }
        }
        if (foundParticle == true) {
          break;
        }
      }
    }
  }

  //-----------------------------------------------------------------------------
  // find the best index in _plottingData to plot (the search that went furthest)
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::findBestPlotIndex() {

    bool helixFound = false;
    int mostHits = 0;
    for (size_t i = 0; i < _plottingData.size(); i++) {
      if (helixFound == true) {
        if (_plottingData[i].becameHelix == false) {
          continue;
        } else {
          if (_plottingData[i].finalNHits > mostHits) {
            mostHits = _plottingData[i].finalNHits;
            _bestPlotIndex = i;
          }
        }
      }
      if (helixFound == false) {
        if (_plottingData[i].becameHelix == true) {
          helixFound = true;
          mostHits = _plottingData[i].finalNHits;
          _bestPlotIndex = i;
        } else {
          if (_plottingData[i].finalNHits > mostHits) {
            mostHits = _plottingData[i].finalNHits;
            _bestPlotIndex = i;
          }
        }
      }
    }
  }

  //-----------------------------------------------------------------------------
  // function to make all the plots in debug mode when runDisplay is on
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::doAllPlots() {

    findBestPlotIndex();
    _xy0->Clear();
    _xy1->Clear();
    _xy2->Clear();
    _phiz2->Clear();
    _bestLineSeg->Clear();
    _bestLineSegResolved->Clear();
    _xy3->Clear();
    _phiz3->Clear();
    _xy4->Clear();
    _phiz4->Clear();
    if (_plottingData[_bestPlotIndex].tcHitsColl.size() > 0) {
      plotXY(0);
    }
    if (_plottingData[_bestPlotIndex].tcHitsColl.size() > 1) {
      plotXY(1);
    }
    if (_plottingData[_bestPlotIndex].tcHitsColl.size() > 2) {
      plotXY(2);
      plotPhiZ(2);
    }
    if (_plottingData[_bestPlotIndex].tcHitsColl.size() > 3) {
      if (_plottingData[_bestPlotIndex].foundLineSegment == true) {
        plotSegment(0);
        plotSegment(1);
      }
      plotXY(3);
      plotPhiZ(3);
    }
    if (_plottingData[_bestPlotIndex].tcHitsColl.size() > 4) {
      plotXY(4);
      plotPhiZ(4);
    }
  }

  //-----------------------------------------------------------------------------
  // function to plot XY view at various stages of logic
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::plotXY(int stage) {

    if (stage == 0) {
      _xy0->cd();
    }
    if (stage == 1) {
      _xy1->cd();
    }
    if (stage == 2) {
      _xy2->cd();
    }
    if (stage == 3) {
      _xy3->cd();
    }
    if (stage == 4) {
      _xy4->cd();
    }

    std::vector<cHit> hits = _plottingData[_bestPlotIndex].tcHitsColl[stage];
    ::LsqSums4 circleFit = _plottingData[_bestPlotIndex].circleFitter[stage];

    // set up MultiGraph that multiple graphs can be added to
    TMultiGraph* mg = new TMultiGraph();
    mg->GetXaxis()->SetLimits(-900, 900);
    mg->GetYaxis()->SetRangeUser(-900, 900);
    mg->GetXaxis()->SetTitle("x");
    mg->GetYaxis()->SetTitle("y");

    // set up graph point for tracker center, helix MC center, and fit center
    float trackerCenterXarray[1] = {0.0};
    float trackerCenterYarray[1] = {0.0};
    TGraph* trackerCenter = new TGraph(1, trackerCenterXarray, trackerCenterYarray);
    trackerCenter->SetMarkerStyle(2);
    trackerCenter->SetMarkerColor(1);
    mg->Add(trackerCenter);
    float mcCenterXarray[1] = {_mcX0};
    float mcCenterYarray[1] = {_mcY0};
    TGraph* mcCenter = new TGraph(1, mcCenterXarray, mcCenterYarray);
    mcCenter->SetMarkerStyle(2);
    mcCenter->SetMarkerColor(4);
    mg->Add(mcCenter);
    mg->Draw("AP");
    double fitCenterXarray[1] = {circleFit.x0()};
    double fitCenterYarray[1] = {circleFit.y0()};
    TGraph* fitCenter = new TGraph(1, fitCenterXarray, fitCenterYarray);
    fitCenter->SetMarkerStyle(2);
    fitCenter->SetMarkerColor(6);
    mg->Add(fitCenter);
    mg->Draw("AP");

    // // draw tracker hollow circle
    TEllipse* hollowCircle = new TEllipse(0.0, 0.0, 380.0, 380.0);
    hollowCircle->SetLineColor(1);
    hollowCircle->SetFillStyle(0);
    hollowCircle->Draw("same");

    // // draw tracker outer edge circle
    TEllipse* outerCircle = new TEllipse(0.0, 0.0, 700.0, 700.0);
    outerCircle->SetLineColor(1);
    outerCircle->SetFillStyle(0);
    outerCircle->Draw("same");

    // draw helix MC circle
    TEllipse* mcCircle = new TEllipse(_mcX0, _mcY0, _mcRadius, _mcRadius);
    mcCircle->SetLineColor(4);
    mcCircle->SetFillStyle(0);
    mcCircle->Draw("same");

    // draw fit circle
    TEllipse* fitCircle =
      new TEllipse(circleFit.x0(), circleFit.y0(), circleFit.radius(), circleFit.radius());
    fitCircle->SetLineColor(6);
    fitCircle->SetFillStyle(0);
    fitCircle->Draw("same");

    // draw calo cluster if it exists, count nParticle points and nBkgPoints
    int nParticlePoints(0), nBkgPoints(0), nParticleSH(0), nBkgSH(0);
    for (size_t i = 0; i < hits.size(); i++) {
      if (stage != 0) {
        if (hits[i].used == false) {
          continue;
        }
      }
      if (hits[i].hitIndice == HitType::STOPPINGTARGET) {
        continue;
      }
      if (hits[i].hitIndice == HitType::CALOCLUSTER) {
        TEllipse* caloCluster =
          new TEllipse(getPos(i).x(), getPos(i).y(), _caloClusterSigma, _caloClusterSigma);
        caloCluster->SetLineColor(807);
        caloCluster->SetFillColorAlpha(807, 0.0);
        caloCluster->Draw("same");
        continue;
      }
      if (hits[i].debugParticle == true) {
        nParticlePoints += 1;
        int hitIndice = hits[i].hitIndice;
        nParticleSH += _chColl->at(hitIndice).nStrawHits();
      }
      if (hits[i].debugParticle == false) {
        nBkgPoints += 1;
        int hitIndice = hits[i].hitIndice;
        nBkgSH += _chColl->at(hitIndice).nStrawHits();
      }
    }

    // set up appropriate number of TLines to draw error bars
    std::vector<TLine*> particleWireLines(nParticlePoints);
    std::vector<TLine*> particlePerpLines(nParticlePoints);
    std::vector<TLine*> bkgWireLines(nBkgPoints);
    std::vector<TLine*> bkgPerpLines(nBkgPoints);

    // add "legend" to plot
    std::vector<std::string> labels = {Form("N_{green} = %d (%d)", nParticlePoints, nParticleSH),
                                       Form("N_{red} = %d (%d)", nBkgPoints, nBkgSH)};
    TLatex* latex = new TLatex();
    latex->SetNDC();
    latex->SetTextFont(43);
    latex->SetTextSize(24);
    latex->SetTextAlign(31);
    float th = 0.05;
    float tx = 0.30;
    float ty = 0.85;
    for (size_t l = 0; l < labels.size(); l++) {
      latex->DrawLatex(tx, ty, labels[l].c_str());
      ty -= th;
    }

    // draw hit error bars along wire and transverse to wire
    int particlePlotIndex = 0;
    int bkgPlotIndex = 0;
    for (size_t i = 0; i < hits.size(); i++) {
      if (stage != 0) {
        if (hits[i].used == false) {
          continue;
        }
      }
      if (hits[i].hitIndice < 0) {
        continue;
      }
      float xPos = getPos(i).x();
      float yPos = getPos(i).y();
      // error bar along wire
      int hitIndice = hits[i].hitIndice;
      float wireDirX = _chColl->at(hitIndice).uDir().x();
      float wireDirY = _chColl->at(hitIndice).uDir().y();
      float wireDirErr = _chColl->at(hitIndice).wireRes();
      float xWire1 = xPos - wireDirX * wireDirErr;
      float xWire2 = xPos + wireDirX * wireDirErr;
      float yWire1 = yPos - wireDirY * wireDirErr;
      float yWire2 = yPos + wireDirY * wireDirErr;
      // error bar transverse to wire
      float perpDirX = _chColl->at(hitIndice).uDir().y();
      float perpDirY = -_chColl->at(hitIndice).uDir().x();
      float perpDirErr = _chColl->at(hitIndice).transRes();
      float xPerp1 = xPos - perpDirX * perpDirErr;
      float xPerp2 = xPos + perpDirX * perpDirErr;
      float yPerp1 = yPos - perpDirY * perpDirErr;
      float yPerp2 = yPos + perpDirY * perpDirErr;
      // initialize lines
      if (hits[i].debugParticle == true) {
        particleWireLines[particlePlotIndex] = new TLine(xWire1, yWire1, xWire2, yWire2);
        particleWireLines[particlePlotIndex]->SetLineColor(8);
        particleWireLines[particlePlotIndex]->Draw("same");
        particlePerpLines[particlePlotIndex] = new TLine(xPerp1, yPerp1, xPerp2, yPerp2);
        particlePerpLines[particlePlotIndex]->SetLineColor(8);
        particlePerpLines[particlePlotIndex]->Draw("same");
        particlePlotIndex++;
      }
      if (hits[i].debugParticle == false) {
        bkgWireLines[bkgPlotIndex] = new TLine(xWire1, yWire1, xWire2, yWire2);
        bkgWireLines[bkgPlotIndex]->SetLineColor(2);
        bkgWireLines[bkgPlotIndex]->Draw("same");
        bkgPerpLines[bkgPlotIndex] = new TLine(xPerp1, yPerp1, xPerp2, yPerp2);
        bkgPerpLines[bkgPlotIndex]->SetLineColor(2);
        bkgPerpLines[bkgPlotIndex]->Draw("same");
        bkgPlotIndex++;
      }
    }
  }

  //-----------------------------------------------------------------------------
  // function to plot phi-z view relative to circle center at various stages of logic
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::plotPhiZ(int stage) {

    if (stage == 2) {
      _phiz2->cd();
    }
    if (stage == 3) {
      _phiz3->cd();
    }
    if (stage == 4) {
      _phiz4->cd();
    }

    std::vector<cHit> hits = _plottingData[_bestPlotIndex].tcHitsColl[stage];
    ::LsqSums2 lineFit = _plottingData[_bestPlotIndex].lineFitter[stage];

    // get values to set range of axis
    float phiMax = 0.0;
    float phiMin = 0.0;
    float zMax = 0.0;
    float zMin = 0.0;
    int hitsConsidered = 0;
    for (size_t i = 0; i < hits.size(); i++) {
      float phi = hits[i].helixPhi + 2 * M_PI * hits[i].helixPhiCorrection;
      float z = getPos(i).z();
      if (hits[i].hitIndice == HitType::STOPPINGTARGET) {
        continue;
      }
      if (hits[i].used == false) {
        continue;
      }
      if (hitsConsidered == 0) {
        phiMax = phi;
        phiMin = phiMax;
        zMax = z;
        zMin = zMax;
        hitsConsidered++;
        continue;
      }
      if (z > zMax) {
        zMax = z;
      }
      if (z < zMin) {
        zMin = z;
      }
      if (phi > phiMax) {
        phiMax = phi;
      }
      if (phi < phiMin) {
        phiMin = phi;
      }
    }

    // set up line for plotting
    float lineSlope = lineFit.dydx();
    float linePhi0 = lineFit.y0();
    float linePhi1 = lineSlope * zMin + linePhi0;
    float linePhi2 = lineSlope * zMax + linePhi0;

    // set up MultiGraph that multiple graphs can be added to
    TMultiGraph* mg = new TMultiGraph();
    mg->GetXaxis()->SetLimits(zMin - 20.0, zMax + 20.0);
    mg->GetYaxis()->SetRangeUser(phiMin - 1.0, phiMax + 1.0);
    mg->GetXaxis()->SetTitle("z");
    mg->GetYaxis()->SetTitle("phi");

    // add dummy point to graph to make sure something is drawn
    float dummyArray1[1] = {0.0};
    float dummyArray2[1] = {0.0};
    TGraph* dummyGraph = new TGraph(1, dummyArray1, dummyArray2);
    dummyGraph->SetMarkerColor(0);
    dummyGraph->SetMarkerSize(0.01);
    mg->Add(dummyGraph);

    // make vectors to make plots for various hit types
    std::vector<float> particlePhiPoints;
    std::vector<float> particleZPoints;
    std::vector<float> particlePhiErrors;
    std::vector<float> particleZErrors;
    std::vector<float> bkgPhiPoints;
    std::vector<float> bkgZPoints;
    std::vector<float> bkgPhiErrors;
    std::vector<float> bkgZErrors;
    std::vector<float> caloPhiPoints;
    std::vector<float> caloZPoints;
    std::vector<float> caloPhiErrors;
    std::vector<float> caloZErrors;
    for (size_t i = 0; i < hits.size(); i++) {
      if (hits[i].used == false) {
        continue;
      }
      if (hits[i].hitIndice == HitType::STOPPINGTARGET) {
        continue;
      }
      float phi = hits[i].helixPhi + 2 * M_PI * hits[i].helixPhiCorrection;
      float helixPhiError = std::sqrt(hits[i].helixPhiError2);
      float z = getPos(i).z();
      float zError = 0.0;
      if (hits[i].hitIndice == HitType::CALOCLUSTER) {
        caloPhiPoints.push_back(phi);
        caloPhiErrors.push_back(helixPhiError);
        caloZPoints.push_back(z);
        caloZErrors.push_back(zError);
        continue;
      }
      if (hits[i].debugParticle == true) {
        particlePhiPoints.push_back(phi);
        particlePhiErrors.push_back(helixPhiError);
        particleZPoints.push_back(z);
        particleZErrors.push_back(zError);
      }
      if (hits[i].debugParticle == false) {
        bkgPhiPoints.push_back(phi);
        bkgPhiErrors.push_back(helixPhiError);
        bkgZPoints.push_back(z);
        bkgZErrors.push_back(zError);
      }
    }

    // make TGraphs
    if (particlePhiPoints.size() > 0) {
      TGraphErrors* particleGraph =
        new TGraphErrors(particlePhiPoints.size(), particleZPoints.data(), particlePhiPoints.data(),
                         particleZErrors.data(), particlePhiErrors.data());
      particleGraph->SetMarkerStyle(20);
      particleGraph->SetMarkerSize(1.1);
      particleGraph->SetMarkerColor(8);
      mg->Add(particleGraph);
    }
    if (bkgPhiPoints.size() > 0) {
      TGraphErrors* bkgGraph =
        new TGraphErrors(bkgPhiPoints.size(), bkgZPoints.data(), bkgPhiPoints.data(),
                         bkgZErrors.data(), bkgPhiErrors.data());
      bkgGraph->SetMarkerStyle(20);
      bkgGraph->SetMarkerSize(1.1);
      bkgGraph->SetMarkerColor(2);
      mg->Add(bkgGraph);
    }
    if (caloPhiPoints.size() > 0) {
      TGraphErrors* caloGraph =
        new TGraphErrors(caloPhiPoints.size(), caloZPoints.data(), caloPhiPoints.data(),
                         caloZErrors.data(), caloPhiErrors.data());
      caloGraph->SetMarkerStyle(20);
      caloGraph->SetMarkerSize(1.1);
      caloGraph->SetMarkerColor(4);
      mg->Add(caloGraph);
    }

    // draw TGraphs
    mg->Draw("AP");

    // add "legend" to plot
    std::vector<std::string> labels = {Form("N_{green} = %d", (int)particlePhiPoints.size()),
                                       Form("N_{red} = %d", (int)bkgPhiPoints.size())};
    TLatex* latex = new TLatex();
    latex->SetNDC();
    latex->SetTextFont(43);
    latex->SetTextSize(24);
    latex->SetTextAlign(31);
    float th = 0.05;
    float tx = 0.30;
    float ty = 0.85;
    for (size_t l = 0; l < labels.size(); l++) {
      latex->DrawLatex(tx, ty, labels[l].c_str());
      ty -= th;
    }

    // add line fit
    if (stage != 2) {
      TLine* phiZFit = new TLine(zMin, linePhi1, zMax, linePhi2);
      phiZFit->SetLineColor(6);
      phiZFit->Draw("same");
    }
  }

  //-----------------------------------------------------------------------------
  // function to plot phi-z segment that leads to best candidate line
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinder::plotSegment(int option) {

    lineInfo line;

    if (option == 0) {
      _bestLineSeg->cd();
      line =
        _plottingData[_bestPlotIndex].lineSegments[_plottingData[_bestPlotIndex].bestLineSegment];
    }
    if (option == 1) {
      _bestLineSegResolved->cd();
      line = _plottingData[_bestPlotIndex]
        .resolvedSegments[_plottingData[_bestPlotIndex].bestLineSegment];
    }

    std::vector<cHit> hits = _plottingData[_bestPlotIndex].tcHitsColl[3];

    // make vectors to make plots for various hit types
    std::vector<float> particlePhiPoints;
    std::vector<float> particleZPoints;
    std::vector<float> particlePhiErrors;
    std::vector<float> particleZErrors;
    std::vector<float> bkgPhiPoints;
    std::vector<float> bkgZPoints;
    std::vector<float> bkgPhiErrors;
    std::vector<float> bkgZErrors;
    for (size_t i = 0; i < line.tcHitsIndices.size(); i++) {
      size_t hitIndex = line.tcHitsIndices[i];
      float phi = hits[hitIndex].helixPhi + 2 * M_PI * line.helixPhiCorrections[i];
      float helixPhiError = std::sqrt(hits[hitIndex].helixPhiError2);
      float z = getPos(hitIndex).z();
      float zError = 0.0;
      if (hits[hitIndex].debugParticle == true) {
        particlePhiPoints.push_back(phi);
        particlePhiErrors.push_back(helixPhiError);
        particleZPoints.push_back(z);
        particleZErrors.push_back(zError);
      }
      if (hits[hitIndex].debugParticle == false) {
        bkgPhiPoints.push_back(phi);
        bkgPhiErrors.push_back(helixPhiError);
        bkgZPoints.push_back(z);
        bkgZErrors.push_back(zError);
      }
    }

    // set up stuff for line
    float lineZMin = line.zMin;
    float lineZMax = line.zMax;
    float lineSlope = line.fitter.dydx();
    float linePhi0 = line.fitter.y0();
    float linePhi1 = lineSlope * lineZMin + linePhi0;
    float linePhi2 = lineSlope * lineZMax + linePhi0;

    // set up MultiGraph that multiple graphs can be added to
    TMultiGraph* mg = new TMultiGraph();
    mg->GetXaxis()->SetLimits(lineZMin - 100.0, lineZMax + 100.0);
    mg->GetYaxis()->SetRangeUser(linePhi1 - 1.0, linePhi2 + 1.0);
    mg->GetXaxis()->SetTitle("z");
    mg->GetYaxis()->SetTitle("phi");

    // add dummy point to graph to make sure something is drawn
    float dummyArray1[1] = {0.0};
    float dummyArray2[1] = {0.0};
    TGraph* dummyGraph = new TGraph(1, dummyArray1, dummyArray2);
    dummyGraph->SetMarkerColor(0);
    dummyGraph->SetMarkerSize(0.01);
    mg->Add(dummyGraph);

    // make TGraphs
    if (particlePhiPoints.size() > 0) {
      TGraphErrors* particleGraph =
        new TGraphErrors(particlePhiPoints.size(), particleZPoints.data(), particlePhiPoints.data(),
                         particleZErrors.data(), particlePhiErrors.data());
      particleGraph->SetMarkerStyle(20);
      particleGraph->SetMarkerSize(1.1);
      particleGraph->SetMarkerColor(8);
      mg->Add(particleGraph);
    }
    if (bkgPhiPoints.size() > 0) {
      TGraphErrors* bkgGraph =
        new TGraphErrors(bkgPhiPoints.size(), bkgZPoints.data(), bkgPhiPoints.data(),
                         bkgZErrors.data(), bkgPhiErrors.data());
      bkgGraph->SetMarkerStyle(20);
      bkgGraph->SetMarkerSize(1.1);
      bkgGraph->SetMarkerColor(2);
      mg->Add(bkgGraph);
    }

    // draw TGraphs
    mg->Draw("AP");

    // add "legend" to plot
    std::vector<std::string> labels = {Form("N_{green} = %d", (int)particlePhiPoints.size()),
                                       Form("N_{red} = %d", (int)bkgPhiPoints.size())};
    TLatex* latex = new TLatex();
    latex->SetNDC();
    latex->SetTextFont(43);
    latex->SetTextSize(24);
    latex->SetTextAlign(31);
    float th = 0.05;
    float tx = 0.30;
    float ty = 0.85;
    for (size_t l = 0; l < labels.size(); l++) {
      latex->DrawLatex(tx, ty, labels[l].c_str());
      ty -= th;
    }

    // add line fit
    TLine* phiZFit = new TLine(lineZMin, linePhi1, lineZMax, linePhi2);
    phiZFit->SetLineColor(6);
    phiZFit->Draw("same");
  }

} // namespace mu2e

using mu2e::AgnosticHelixFinder;
DEFINE_ART_MODULE(AgnosticHelixFinder)
