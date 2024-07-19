//
// Object to perform helix fit to straw hits
//
//
#ifndef CalHelixFinderAlg_HH
#define CalHelixFinderAlg_HH

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

// data
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"

// BaBar
#include "BTrk/TrkBase/TrkErrCode.hh"
//root
//#include "TString.h"

#include "Offline/Mu2eUtilities/inc/LsqSums2.hh"
#include "Offline/Mu2eUtilities/inc/LsqSums4.hh"
// #include "CalPatRec/inc/CalTimePeak.hh"
//#include "CalPatRec/inc/CalHelixPoint.hh"
#include "Offline/CalPatRec/inc/CalHelixFinderData.hh"

class TH1F;

namespace mu2e {
  class Calorimeter;
  class Tracker;
//-----------------------------------------------------------------------------
// output struct
//-----------------------------------------------------------------------------
  class CalHelixFinderAlg {
  public:
    enum { kMaxNHits = 10000 } ;

    struct Config {
      fhicl::Atom<int>            diag{                fhicl::Name("diagLevel"),                    fhicl::Comment("Diag"),0 };
      fhicl::Atom<int>            debug{               fhicl::Name("debugLevel"),                   fhicl::Comment("Debug"),0 };
      fhicl::Atom<int>            debug2{              fhicl::Name("debugLevel2"),                  fhicl::Comment("Debug2"),0 };
      fhicl::Sequence<std::string> hsel{               fhicl::Name("HelixFitSelectionBits"),                    fhicl::Comment("Good Hit Selection") };
      fhicl::Sequence<std::string> bkgsel{             fhicl::Name("BackgroundSelectionBits"),                  fhicl::Comment("Background Hit Selection") };
      fhicl::Atom<float>          maxHitEnergy{        fhicl::Name("maxElectronHitEnergy"),            fhicl::Comment("MaxElectronHitEnergy") };
      fhicl::Atom<int>            minNHits{            fhicl::Name("minNHit"),                fhicl::Comment("Min NHits") };
      fhicl::Atom<float>          absMpDfDz{           fhicl::Name("mostProbableDfDz"),               fhicl::Comment("Most Probable DfDz") };
      fhicl::Atom<int>            initDfDz{            fhicl::Name("initDfDz"),                fhicl::Comment("Initial DfDz") };
      fhicl::Atom<float>            dzOverHelPitchCut{   fhicl::Name("dzOverHelPitchCut"),       fhicl::Comment("Cut on Ratio Between Dz and HelPitch") };
      fhicl::Atom<float>          maxDfDz{             fhicl::Name("maxDfDz"),                 fhicl::Comment("Max DfDz") };
      fhicl::Atom<float>          minDfDz{             fhicl::Name("minDfDz"),                 fhicl::Comment("Min DfDz") };
      fhicl::Atom<float>          sigmaPhi{            fhicl::Name("sigmaPhi"),                fhicl::Comment("Sigma Phi") };
      fhicl::Atom<float>          weightXY{            fhicl::Name("weightXY"),                fhicl::Comment("Weight XY") };
      fhicl::Atom<int>            targetcon{           fhicl::Name("targetconsistent"),               fhicl::Comment("Target Consistent") };
      fhicl::Atom<float>          weightZPhi{          fhicl::Name("weightZPhi"),              fhicl::Comment("Weight ZPhi") };
      fhicl::Atom<float>          weight3D{            fhicl::Name("weight3D"),                fhicl::Comment("Weight 3D") };
      fhicl::Atom<float>          maxXDPhi{            fhicl::Name("maxXDPhi"),                fhicl::Comment("MaxXDPhi") };
      fhicl::Atom<float>          maxPanelToHelixDPhi{ fhicl::Name("maxPanelToHelixDPhi"),     fhicl::Comment("Max Panel to Helix DPhi") };
      fhicl::Atom<float>          distPatRec{          fhicl::Name("distPatRec"),              fhicl::Comment("Dist Pat Rec") };
      fhicl::Atom<float>          minDeltaNShPatRec{   fhicl::Name("minDeltaNShPatRec"),       fhicl::Comment("Min Delta NSh Pat Rec") };
      fhicl::Atom<float>          mindist{             fhicl::Name("mindist"),                 fhicl::Comment("Min Dist") };
      fhicl::Atom<float>          pmin{                fhicl::Name("minP"),                    fhicl::Comment("P Min") };
      fhicl::Atom<float>          pmax{                fhicl::Name("maxP"),                    fhicl::Comment("P Max") };
      fhicl::Atom<float>          tdmin{               fhicl::Name("minAbsTanDip"),                   fhicl::Comment("Min Abs TanDip") };
      fhicl::Atom<float>          tdmax{               fhicl::Name("maxAbsTanDip"),                   fhicl::Comment("Max Abs TanDip") };
      fhicl::Atom<bool>           xyweights{           fhicl::Name("xyWeights"),               fhicl::Comment("XY Weights") };
      fhicl::Atom<bool>           zweights{            fhicl::Name("zWeights"),                fhicl::Comment("Z Weights") };
      fhicl::Atom<bool>           filter{              fhicl::Name("filter"),                  fhicl::Comment("Filter") };
      fhicl::Atom<bool>           plotall{             fhicl::Name("plotall"),                 fhicl::Comment("Plot All") };
      fhicl::Atom<bool>           usetarget{           fhicl::Name("usetarget"),               fhicl::Comment("Use Target") };
      fhicl::Atom<float>          maxZTripletSearch{   fhicl::Name("maxZTripletSearch"),       fhicl::Comment("Max Z Triplet Search") };
      fhicl::Atom<float>          bz{                  fhicl::Name("bz"),                      fhicl::Comment("Value of field z component")};
      fhicl::Atom<int>            nHitsMaxPerPanel{    fhicl::Name("nHitsMaxPerPanel"),        fhicl::Comment("Max NHits Per Panel") };
      fhicl::Atom<float>          hitChi2Max{          fhicl::Name("hitChi2Max"),              fhicl::Comment("Hit Chi2 Max") };
      fhicl::Atom<float>          chi2xyMax{           fhicl::Name("chi2xyMax"),               fhicl::Comment("Chi2XY Max") };
      fhicl::Atom<float>          chi2zphiMax{         fhicl::Name("chi2zphiMax"),             fhicl::Comment("Chi2ZPhi Max") };
      fhicl::Atom<float>          chi2hel3DMax{        fhicl::Name("chi2hel3DMax"),            fhicl::Comment("Chi2 Hel3D Max") };
      fhicl::Atom<float>          dfdzErr{             fhicl::Name("dfdzErr"),                 fhicl::Comment("DfDz Error") };
      fhicl::Atom<float>          maxNHitsRatio{       fhicl::Name("maxNHitsRatio"),           fhicl::Comment("Max NHits Ratio") };
      fhicl::Atom<float>          minArea{             fhicl::Name("minArea"),                 fhicl::Comment("Minimum triplet area") };
    };
//-----------------------------------------------------------------------------
// data members
//-----------------------------------------------------------------------------
    const Tracker*            _tracker;
    const Calorimeter*         _calorimeter;

    //    const CalTimePeak*         fTimePeak;
    //    const TimeCluster*         fTimeCluster; //needed for debugging

    float                     fCaloTime;
    float                     fCaloX;
    float                     fCaloY;
    float                     fCaloZ;

    //    std::vector<CalHelixPoint> _xyzp;        // normally includes only hits from the time peak
//-----------------------------------------------------------------------------
// for diagnostics purposes save several states of _xyzp (only if _diag > 0)
//
// [0]: copy of initial _xyzp state
// [1]: after filterDist
// [2]: doPatternRecognition: after rescueHitsBeforeSeed
// [3]: doPatternRecognition: after rescueHits  (inside doPatternRecognition)
// [4]: doPatternRecognition: after final findZ (inside doPatternRecognition)
// [5]: doPatternRecognition: after doPatternRecognition
//
// each XYZPHack point has an 'outlier' flag telling if it belongs to the track
//-----------------------------------------------------------------------------
    struct SaveResults_t {
      //      std::vector<CalHelixPoint> _xyzp;    // points with used flags
      CalHelixFinderData         _helix;   // helix, found at this stage
    };

    struct PhiZFitInfo_t     {
      HitInfo_t          seedIndex;
      float             dfdz;
      float             phi0;
      int                faceId;
      float             weight;
      int                useInteligentWeight;
      float             dz;
      float             zlast;
      int                nPointsRemoved;
      int                doCleanUp;
      const char*        banner;
    };

    SaveResults_t        _results[6];      // diagnostic buffers

    int                  fUseDefaultDfDz;

    int                  _diag;
    int                  _debug;
    int                  _debug2;
    // int                  _smartTag;     //flag used to test addiotional layer of rejection after the search for the "best triplet"
    StrawHitFlag         _hsel;         // good hit selection
    StrawHitFlag         _bkgsel;       // background hit selection
    float               _maxHitEnergy; //
    int                  _minNHits;     // minimum # of hits for a helix candidate
                                        // 2014-03-10 Gianipez and P. Murat: limit
                                        // the dfdz value in the pattern-recognition stage
    float               _mpDfDz;
    float               _absMpDfDz;         // absolute value of most probable expected dphi/dz
    int                 _initDfDz;
    float               _dzOverHelPitchCut; //cut on the ratio between the Dz and the predicted helix-pitch used in ::findDfDz(...)
    float               _maxDfDz;
    float               _minDfDz;
    float               _sigmaPhi;     // hit phi resolution (wrt the trajectory axis, assume R=25-30cm)
    float               _weightXY;     // scale factor for makeing the xy-chi2 with a mean close to 1
    int                  _targetcon;    //require or not the circle fit to inlcude the Stopping Target center
    float               _weightZPhi;
    float               _weight3D;
    float               _maxXDPhi;     // max normalized hit residual in phi (findRZ)
    float               _maxPanelToHelixDPhi;  // max dphi between the helix prediction and a given tracker plane

    float               _hdfdz;        // estimated d(phi)/dz value
    float               _sdfdz;        // estimated d(phi)/dz error
    float               _hphi0;
                                        // 201-03-31 Gianipez added for changing the value of the
                                        // squared distance requed bewtween a straw hit and its predicted
                                        // position used in the patter recognition procedure
    float               _distPatRec;
    int                 _minDeltaNShPatRec;  //minimum number of additional StrawHits required in
                                              //the findTrack function to set the new Helix

    float               _mindist;       // minimum distance between points used in circle initialization
    float               _pmin, _pmax;   // range of total momentum
    float               _tdmin, _tdmax; // range of abs(tan(dip)
    bool                 _xyweights;
    bool                 _zweights;      // weight points by estimated errors
    bool                 _filter;        // filter hits
    bool                 _plotall;       // plot also failed fits
    bool                 _usetarget;     // constrain to target when initializing
    float                _maxZTripletSearch; //maximum z allowed for the hit used to search the best triplet
    mutable float       _bz;            // cached value of Field Z component at the tracker origin
//-----------------------------------------------------------------------------
// cached value of radius and pitch sign: these depend on the particle type
// and direction
//-----------------------------------------------------------------------------
    float    _rmin, _rmax, _smin, _smax, _dfdzsign;
//-----------------------------------------------------------------------------//
// store the paramters value of the most reliable track candidate
//-----------------------------------------------------------------------------//
    int       _nHitsMaxPerPanel;

                                         // thresholds for the worst hit chi2, total XY and ZPhi fit chi2's
    float    _hitChi2Max;
    float    _chi2xyMax;
    float    _chi2zphiMax;
    float    _chi2hel3DMax;

    // indices, distance from prediction and distance along z axis from the seeding hit
    // of the hits found in the pattern recognition

    int       _phiCorrectedDefined;

    float    _dfdzErr;                 // error on dfdz by ::findDfDz
    float    _maxNHitsRatio;
    float    _minarea2;
//-----------------------------------------------------------------------------
// checkpoints, used for debugging
//-----------------------------------------------------------------------------
    int       _findTrackLoopIndex;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  public:
                                        // parameter set should be passed in on construction
    //explicit CalHelixFinderAlg(fhicl::ParameterSet const&);
    explicit CalHelixFinderAlg(const Config& config);
    virtual ~CalHelixFinderAlg();
                                        // cached bfield accessor
    float bz() const;

    void   calculateDfDz    (float phi0, float phi1, float z0,  float z1, float &dfdz);
    void   calculateDphiDz_2(CalHelixFinderData& Helix,HitInfo_t HitIndex, int NHits, float X0, float Y0, float& DphiDz);

    //projects the straw hit error along the radial direction of the circle-helix
    float calculateWeight     (const mu2e::ComboHit&  Hit,
                                // const CLHEP::Hep3Vector& HitPos,
                                // const CLHEP::Hep3Vector& StrawDir,
                                const XYZVectorF& HelCenter,
                                float                   Radius);

    float calculatePhiWeight  (const ComboHit&  Hit,
                                // const XYZVectorF& HitPos   ,
                                // const XYZVectorF& StrawDir ,
                                const XYZVectorF& HelCenter,
                                float                   Radius   ,
                                int                      Print    ,
                                const char*              Banner=NULL);

    void   resolve2PiAmbiguity (ComboHit* Hit, XYZVectorF& Center, float &Phi_ref, float &DPhi);

    //calculates the residual along the radial direction of the helix-circle
    float calculateRadialDist (const XYZVectorF& HitPos,
                                const XYZVectorF& HelCenter,
                                float                   Radius);

    bool   calculateTrackParameters(const XYZVectorF& p1,
                                    const XYZVectorF& p2,
                                    const XYZVectorF& p3,
                                    XYZVectorF&       Center,
                                    float&       Radius,
                                    float&       Phi0,
                                    float&       TanLambda);

    void   countUsedHits           (CalHelixFinderData& Helix,
                                    HitInfo_t           SeedIndex,
                                    int&                NComboHits,
                                    int&                NPoints);
    static float deltaPhi         (float phi1, float phi2);

   // returns the index of the hit which provides the highest contribute to the chi2
    void   cleanUpWeightedCircleFit(CalHelixFinderData& Helix,
                                    HitInfo_t           SeedIndex,
                                    HitInfo_t&          Iworst);

    bool   doLinearFitPhiZ          (CalHelixFinderData& Helix,
                                     HitInfo_t           SeedIndex,
                                     int                 UseInteligentWeight=0,
                                     int                 DoCleanUp =1);

    void   addCaloClusterToFitPhiZ  (CalHelixFinderData& Helix);

    void   findGoodFaceHitInFitPhiZ (CalHelixFinderData& Helix,
                                     PhiZFitInfo_t&      phiZInfo,
                                     HitInfo_t&          GoodHit,
                                     float&              HitChi2);

    void   findWorstChi2HitInFitPhiZ(CalHelixFinderData& Helix,
                                     PhiZFitInfo_t&      phiZInfo,
                                     HitInfo_t&          WorstHit,
                                     float&             HitChi2);


    void   findWorstResidHitInFitPhiZ(CalHelixFinderData& Helix,
                                      PhiZFitInfo_t&      phiZInfo,
                                      HitInfo_t&          WorstHit,
                                      float&             HitChi2);

    // void   doCleanUpFitPhiZ         (CalHelixFinderData& Helix,
    //                                      );
   //perfoms the weighted circle fit, update the helix parameters (HelicCenter, Radius) and
    // fills the vector Weights which holds the calculated weights of the hits
    void   doWeightedCircleFit     (CalHelixFinderData& Helix,
                                    HitInfo_t           SeedIndex,
                                    XYZVectorF&             HelCenter,
                                    float&             Radius,
                                    int                 Print=0,
                                    const char*         Banner=NULL);

    void  doPatternRecognition(CalHelixFinderData& mytrk);


    //performs the search of the best triplet
    void  searchBestTriplet   (CalHelixFinderData& Helix, CalHelixFinderData& TmpHelix, int UseMPVdfdz=0);

    void  defineHelixParams   (CalHelixFinderData& Helix) const;

    // TH1F* hDist() {return _hDist;}

    int   isHitUsed(int index);

    void  fillFaceOrderedHits          (CalHelixFinderData& Helix);
    // void filterDist                   (CalHelixFinderData& Helix);
    void  filterUsingPatternRecognition(CalHelixFinderData& Helix);
    void  setCaloCluster               (CalHelixFinderData& Helix);
    bool  fitHelix                     (CalHelixFinderData& Helix);
    // bool findHelix                    (CalHelixFinderData& Helix, const CalTimePeak* TimePeak);
    bool findHelix                    (CalHelixFinderData& Helix);
    int  findDfDz                     (CalHelixFinderData& Helix, HitInfo_t SeedIndex, int  Diag_flag=0);
    int  findDfDz_1                   (CalHelixFinderData& Helix, HitInfo_t SeedIndex, int  Diag_flag=0);
    int  findDfDz_2                   (CalHelixFinderData& Helix, HitInfo_t SeedIndex, int  Diag_flag=0);
    void findTrack                    (HitInfo_t&         SeedIndex,
                                       CalHelixFinderData& Helix,
                                       int                 UseMPVdfdz     = 0);

    // float ApproxAtan                  (float z);
    // float polyAtan2                   (float y, float x);
    bool isFaceUsed                   (CalHelixFinderData& Helix, FaceZ_t* face);
//-----------------------------------------------------------------------------
// setters
//-----------------------------------------------------------------------------
    void  setTracker    (const Tracker*    Tracker) { _tracker     = Tracker; }
    void  setCalorimeter(const Calorimeter* Cal    ) { _calorimeter = Cal    ; }
//-----------------------------------------------------------------------------
// diagnostics
//-----------------------------------------------------------------------------
    void   printInfo            (CalHelixFinderData& Helix);
    void   printXYZP            (CalHelixFinderData& Helix);

    int    refineHelixParameters(CalHelixFinderData& Helix,
                                 HitInfo_t          SeedIndex,
                                 const char*         Banner=NULL,
                                 int                 Print=0);

                                        // 12-10-2013 Gianipez: new pattern recognition functions
    void   rescueHitsBeforeSeed (CalHelixFinderData&  Helix);

    void   rescueHits           (CalHelixFinderData&  Helix, HitInfo_t SeedIndex   ,
                                 int UsePhiResiduals = 0);

    // void   resolve2PiAmbiguity  (CalHelixFinderData& Helix,const XYZVectorF& Center, float DfDz, float Phi0);

    void   resetTrackParamters  ();
//-----------------------------------------------------------------------------
// save intermediate results in diagnostics mode
//-----------------------------------------------------------------------------
    void   saveResults                    (CalHelixFinderData&         Helix, int Index);

    // void   fillHelixDiag                  (CalHelixFinderData&         Helix);

    void   searchWorstHitWeightedCircleFit(CalHelixFinderData& Helix,
                                           HitInfo_t          SeedIndex,
                                           // int*               IdVec,
                                           const XYZVectorF& HelCenter,
                                           float&             Radius,
                                           // float*            Weights,
                                           HitInfo_t&         Iworst ,
                                           float&             HitChi2Worst);



  };
}
#endif
