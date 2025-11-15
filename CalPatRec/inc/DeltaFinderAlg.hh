//
// Compton electron finder algorithm
//
//
#ifndef CalPatRec_DeltaFinderAlg_hh
#define CalPatRec_DeltaFinderAlg_hh

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/ParameterSet.h"

#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"

#include "Offline/Mu2eUtilities/inc/LsqSums2.hh"
#include "Offline/Mu2eUtilities/inc/StopWatch.hh"

#include "Offline/CalPatRec/inc/DeltaFinder_types.hh"

namespace fhicl {
  class ParameterSet;
}

namespace mu2e {
  class Calorimeter;
  class Tracker;

  using namespace DeltaFinderTypes;

  class DeltaFinderAlg {
  public:
//-----------------------------------------------------------------------------
// algorithm talk-to parameters
//-----------------------------------------------------------------------------
    struct Config {
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<int>             debugLevel        {Name("debugLevel"        ), Comment("debug level"                 ) };
      fhicl::Atom<int>             diagLevel         {Name("diagLevel"         ), Comment("diag level"                  ) };
      fhicl::Atom<int>             printErrors       {Name("printErrors"       ), Comment("print errors"                ) };
      fhicl::Atom<int>             flagProtonHits    {Name("flagProtonHits"    ), Comment("1: flag proton hits"         ) };
      fhicl::Atom<int>             mergePC           {Name("mergePC"           ), Comment("1: merge proton candidates"  ) };
      fhicl::Atom<int>             pickupProtonHits  {Name("pickupProtonHits"  ), Comment("1: pickup single proton hits") };
      fhicl::Atom<float>           timeBin           {Name("timeBin"           ), Comment("time bin for storing hits"   ) };
      fhicl::Atom<float>           maxDeltaEDep      {Name("maxDeltaEDep"      ), Comment("max delta candidate  eDep"   ) };
      fhicl::Atom<float>           maxSeedEDep       {Name("maxSeedEDep"       ), Comment("max seed eDep"               ) };
      fhicl::Atom<float>           minProtonSeedEDep {Name("minProtonSeedEDep" ), Comment("min proton seed eDep"        ) };
      fhicl::Atom<float>           minProtonHitEDep  {Name("minProtonHitEDep"  ), Comment("min proton hit eDep"         ) };
      fhicl::Atom<int>             minNSeeds         {Name("minNSeeds"         ), Comment("min N seeds in a delta cand" ) };
      fhicl::Atom<int>             minDeltaNHits     {Name("minDeltaNHits"     ), Comment("min N combo  hits in a delta") };
      fhicl::Atom<float>           maxEleHitEnergy   {Name("maxEleHitEnergy"   ), Comment("max electron hit energy"     ) };
      fhicl::Atom<float>           maximumTime       {Name("maximumTime"       ), Comment("maximum time"                ) };
      fhicl::Atom<float>           maxHitSeedDt      {Name("maxHitSeedDt"      ), Comment("max DT(hit-seed)"            ) };
      fhicl::Atom<float>           maxChi2Seed       {Name("maxChi2Seed"       ), Comment("max seed chi2 (stereo)"      ) };
      // fhicl::Atom<float>           scaleTwo          {Name("scaleTwo"          ), Comment("scale factor chi2 (stereo)"  ) };
      fhicl::Atom<float>           maxChi2Par        {Name("maxChi2Par"        ), Comment("max chi2 (parallel)"         ) };
      fhicl::Atom<float>           maxChi2Perp       {Name("maxChi2Perp"       ), Comment("max chi2 (perp to the wire)" ) };
      // fhicl::Atom<float>           maxChi2All        {Name("maxChi2All"        ), Comment("max chi2 (all)"              ) };
      fhicl::Atom<float>           maxChi2SeedDelta  {Name("maxChi2SeedDelta"  ), Comment("max chi2 (seed-delta)"       ) };
      fhicl::Atom<float>           rCore             {Name("rCore"             ), Comment("core radius"                 ) };
      fhicl::Atom<int>             maxGap            {Name("maxGap"            ), Comment("max Gap"                     ) };
      fhicl::Atom<float>           sigmaR            {Name("sigmaR"            ), Comment("sigmaR"                      ) };
      fhicl::Atom<float>           maxDriftTime      {Name("maxDriftTime"      ), Comment("maxDriftTime"                ) };
      fhicl::Atom<float>           maxSeedDt         {Name("maxSeedDt"         ), Comment("maxSeedDt"                   ) };
      fhicl::Atom<float>           maxHitDt          {Name("maxHitDt"          ), Comment("maxHitDt"                    ) };
      // fhicl::Atom<float>           maxStrawDt        {Name("maxStrawDt"        ), Comment("max straw Dt"                ) };
      fhicl::Atom<float>           maxDtDs           {Name("maxDtDs"           ), Comment("max Dt/Dstation"             ) };
      fhicl::Atom<float>           maxDtDc           {Name("maxDtDc"           ), Comment("max deltaT between deltas"   ) };
      fhicl::Atom<int>             testOrder         {Name("testOrder"         ), Comment("if 1, test order"            ) };
      fhicl::Atom<bool>            testHitMask       {Name("testHitMask"       ), Comment("if true, test hit mask"      ) };
      fhicl::Sequence<std::string> goodHitMask       {Name("goodHitMask"       ), Comment("good hit mask"               ) };
      fhicl::Sequence<std::string> bkgHitMask        {Name("bkgHitMask"        ), Comment("background hit mask"         ) };
    };
  public:
//-----------------------------------------------------------------------------
// data members
//-----------------------------------------------------------------------------
    Data_t*         _data;

    int             _flagProtonHits;       // if 1, find protons and flag proton hits
    int             _mergePC;              //
    int             _pickupProtonHits;     //
    float           _timeBin;              // binning
    // float           _minHitTime;           // min hit time
    float           _maxDeltaEDep;         //
    float           _maxSeedEDep;          //
    float           _minProtonSeedEDep;    //
    float           _minProtonHitEDep;     //
    int             _minNSeeds;            // min number of seeds in the delta electron cluster
    int             _minDeltaNHits;        // min number of hits of a delta candidate
    float           _maxEleHitEnergy;      //
    // float           _minT;
    float           _maxT;
    float           _maxHitSeedDt;         //
    float           _maxChi2Seed;          //
    // float           _scaleTwo;
    // float           _maxChi2Neighbor;      //
    float           _maxChi2Par;           //
    float           _maxChi2Perp;          //
    // float           _maxChi2All;           // max chi2/N of a seed
    float           _maxChi2SeedDelta;     // max delta-seed chi2 for adding a seed to a delta
    float           _rCore;                // represents "core radius" of the conversion trajectory
    int             _maxGap;
    float           _sigmaR;
    float           _sigmaR2;              // _sigmaR^2
    float           _maxDriftTime;
    float           _maxSeedDt;            // +/- SeedDt is the time window for checking duplicate seeds
    float           _maxHitDt;
    // float           _maxStrawDt;
    float           _maxDtDs;              // low-P electron travel time between two stations
    float           _maxDtDc;              // max deltaT between two delta candiates

    // int             _writeComboHits;       // write (filtered ?) combo hits
    // int             _writeStrawHitFlags;

    int             _debugLevel;
    int             _diagLevel;
    int             _printErrors;
    int             _testOrder;
    int             _doTiming;

    bool            _testHitMask;
    StrawHitFlag    _goodHitMask;
    StrawHitFlag    _bkgHitMask;

    std::shared_ptr<StopWatch> _watch;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
    DeltaFinderAlg() {}

    explicit     DeltaFinderAlg(const fhicl::Table<DeltaFinderAlg::Config>& config, Data_t* _data);

    float        flagProtonHits()    { return _flagProtonHits; }
    float        printErrors   ()    { return _printErrors   ; }
    float        timeBin       ()    { return _timeBin       ; }
//-----------------------------------------------------------------------------
// other functions
//-----------------------------------------------------------------------------
    int          checkDuplicates     (int Station,
                                      int Face1, const HitData_t* Hit1,
                                      int Face2, const HitData_t* Hit2);

    void         completeSeed        (DeltaSeed* Seed);

    void         findSeeds           (int Station, int Face);
    void         findSeeds           ();
    void         linkDeltaSeeds      ();                        // do it in upstream direction
    int          mergeDeltaCandidates();

    int          orderHits           ();
    void         pruneSeeds          (int Station);
    int          recoverMissingHits  ();
    int          recoverSeed         (DeltaCandidate* Delta, int LastStation, int Station);
    int          recoverStation      (DeltaCandidate* Delta, int LastStation, int Station, int UseUsedHits, int RecoverSeeds);
    void         run                 ();
//-----------------------------------------------------------------------------
// chi2 calculation. Split chi^2 into parallel and perpendicular to the wire components
//-----------------------------------------------------------------------------
    void         deltaChi2(DeltaCandidate* Delta, float Xc, float Yc, float& Chi2Par, float& Chi2Perp);
    void         seedChi2 (DeltaSeed*      Seed , float Xc, float Yc, float& Chi2Par, float& Chi2Perp);
//-----------------------------------------------------------------------------
// make it inline
//-----------------------------------------------------------------------------
    void         hitChi2  (const HitData_t* Hd  , float Xc, float Yc, float& Chi2Par, float& Chi2Perp) {
      float dx        = Hd->fX-Xc;
      float dy        = Hd->fY-Yc;

      float dxy_dot_w = dx*Hd->fWx+dy*Hd->fWy;
      Chi2Par         = (dxy_dot_w*dxy_dot_w)/(_sigmaR2+Hd->fSigW2);

      float dxy_dot_n = dx*Hd->fWy-dy*Hd->fWx;
      float drr       = fmax(fabs(dxy_dot_n)-_rCore,0);
      Chi2Perp        = (drr*drr)/_sigmaR2;
    }
//-----------------------------------------------------------------------------
// proton section of the algorithm
//-----------------------------------------------------------------------------
    int          createProtonCandidates       ();
    int          findProtons                  ();
    int          mergeNonOverlappingCandidates(std::vector<ProtonCandidate*>* Pc);
    int          mergeProtonCandidates        ();
    int          prepareProtonHits            ();
    int          recoverMissingProtonHits     ();
    int          resolveProtonOverlaps        (std::vector<ProtonCandidate*>* Pc);
  };
}
#endif
