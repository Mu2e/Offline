//
// Compton electron finder algorithm
//
//
#ifndef __DeltaFinderAlg_hh__
#define __DeltaFinderAlg_hh__

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/ParameterSet.h"

#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"

#include "Offline/Mu2eUtilities/inc/LsqSums2.hh"

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
      fhicl::Atom<float>           timeBin           {Name("timeBin"           ), Comment("time bin for storing hits"   ) };
      fhicl::Atom<float>           minHitTime        {Name("minHitTime"        ), Comment("min hit time"                ) };
      fhicl::Atom<float>           maxDeltaEDep      {Name("maxDeltaEDep"      ), Comment("max delta candidate  eDep"   ) };
      fhicl::Atom<float>           maxSeedEDep       {Name("maxSeedEDep"       ), Comment("max seed eDep"               ) };
      fhicl::Atom<float>           minProtonSeedEDep {Name("minProtonSeedEDep" ), Comment("min proton seed eDep"        ) };
      fhicl::Atom<int>             minNSeeds         {Name("minNSeeds"         ), Comment("min N seeds in a delta cand" ) };
      fhicl::Atom<int>             minDeltaNHits     {Name("minDeltaNHits"     ), Comment("min N combo  hits in a delta") };
      fhicl::Atom<float>           maxEleHitEnergy   {Name("maxEleHitEnergy"   ), Comment("max electron hit energy"     ) };
      fhicl::Atom<float>           minimumTime       {Name("minimumTime"       ), Comment("minimum time"                ) };
      fhicl::Atom<float>           maximumTime       {Name("maximumTime"       ), Comment("maximum time"                ) };
      fhicl::Atom<float>           maxHitSeedDt      {Name("maxHitSeedDt"      ), Comment("max DT(hit-seed)"            ) };
      fhicl::Atom<float>           maxChi2Seed       {Name("maxChi2Seed"       ), Comment("max seed chi2 (stereo)"      ) };
      fhicl::Atom<float>           scaleTwo          {Name("scaleTwo"          ), Comment("scale factor chi2 (stereo)"  ) };
      fhicl::Atom<float>           maxChi2Radial     {Name("maxChi2Radial"     ), Comment("max chi2 (radial)"           ) };
      fhicl::Atom<float>           maxChi2All        {Name("maxChi2All"        ), Comment("max chi2 (all)"              ) };
      fhicl::Atom<float>           maxChi2SeedDelta  {Name("maxChi2SeedDelta"  ), Comment("max chi2 (seed-delta)"       ) };
      fhicl::Atom<float>           rCore             {Name("rCore"             ), Comment("core radius"                 ) };
      fhicl::Atom<int>             maxGap            {Name("maxGap"            ), Comment("max Gap"                     ) };
      fhicl::Atom<float>           sigmaR            {Name("sigmaR"            ), Comment("sigmaR"                      ) };
      fhicl::Atom<float>           maxDriftTime      {Name("maxDriftTime"      ), Comment("maxDriftTime"                ) };
      fhicl::Atom<float>           maxSeedDt         {Name("maxSeedDt"         ), Comment("maxSeedDt"                   ) };
      fhicl::Atom<float>           maxHitDt          {Name("maxHitDt"          ), Comment("maxHitDt"                    ) };
      fhicl::Atom<float>           maxStrawDt        {Name("maxStrawDt"        ), Comment("max straw Dt"                ) };
      fhicl::Atom<float>           maxDtDs           {Name("maxDtDs"           ), Comment("max Dt/Dstation"             ) };
      fhicl::Atom<float>           maxDtDc           {Name("maxDtDc"           ), Comment("max deltaT between deltas"   ) };
      fhicl::Atom<int>             testOrder         {Name("testOrder"         ), Comment("if 1, test order"            ) };
      fhicl::Atom<bool>            testHitMask       {Name("testHitMask"       ), Comment("if true, test hit mask"      ) };
      fhicl::Sequence<std::string> goodHitMask       {Name("goodHitMask"       ), Comment("good hit mask"               ) };
      fhicl::Sequence<std::string> bkgHitMask        {Name("bkgHitMask"        ), Comment("background hit mask"         ) };
      fhicl::Atom<int>             updateSeedCOG     {Name("updateSeedCOG"     ), Comment("if 1, update seed COG"       ) };
    };
  public:
//-----------------------------------------------------------------------------
// data members
//-----------------------------------------------------------------------------
    Data_t*         _data;

    float           _timeBin;              // binning
    float           _minHitTime;           // min hit time
    float           _maxDeltaEDep;         //
    float           _maxSeedEDep;          //
    float           _minProtonSeedEDep;    //
    int             _minNSeeds;            // min number of seeds in the delta electron cluster
    int             _minDeltaNHits;        // min number of hits of a delta candidate
    float           _maxEleHitEnergy;      //
    float           _minT;
    float           _maxT;
    float           _maxHitSeedDt;         //
    float           _maxChi2Seed;          //
    float           _scaleTwo;
    float           _maxChi2Neighbor;      //
    float           _maxChi2Radial;        //
    float           _maxChi2All;           // max chi2/N of a seed
    float           _maxChi2SeedDelta;     // max delta-seed chi2 for adding a seed to a delta
    float           _rCore;                // represents "core radius" of the conversion trajectory
    int             _maxGap;
    float           _sigmaR;
    float           _sigmaR2;              // _sigmaR^2
    float           _maxDriftTime;
    float           _maxSeedDt;            // +/- SeedDt is the time window for checking duplicate seeds
    float           _maxHitDt;
    float           _maxStrawDt;
    float           _maxDtDs;              // low-P electron travel time between two stations
    float           _maxDtDc;              // max deltaT between two delta candiates

    int             _writeComboHits;       // write (filtered ?) combo hits
    int             _writeStrawHitFlags;

    int             _debugLevel;
    int             _diagLevel;
    int             _printErrors;
    int             _testOrder;

    bool            _testHitMask;
    StrawHitFlag    _goodHitMask;
    StrawHitFlag    _bkgHitMask;

    int             _updateSeedCOG;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
    DeltaFinderAlg() {}

    explicit     DeltaFinderAlg(const fhicl::Table<DeltaFinderAlg::Config>& config, Data_t* _data);

    int          checkDuplicates     (int Station,
                                      int Face1, const HitData_t* Hit1,
                                      int Face2, const HitData_t* Hit2);

    void         completeSeed        (DeltaSeed* Seed);
    void         connectSeeds        ();                        // do it in upstream direction
    void         findSeeds           (int Station, int Face);
    void         findSeeds           ();
    int          mergeDeltaCandidates();

    int          orderHits           ();
//-----------------------------------------------------------------------------
// custom comparator (less) to sort hits in ascending order
// order hits in time, not correctedTime: TOT has tails...
//-----------------------------------------------------------------------------
    // static bool  less             (const ComboHit*& a, const ComboHit*& b) {
    //   return a->time() < b->time();
    // }

    void         pruneSeeds          (int Station);
    int          recoverMissingHits  ();
    int          recoverSeed         (DeltaCandidate* Delta, int LastStation, int Station);
    int          recoverStation      (DeltaCandidate* Delta, int LastStation, int Station, int UseUsedHits, int RecoverSeeds);
    void         run                 ();
                                        // returns chi2^2 sums, not normalized to the number of hits.
                                        // useful when adding a hit to the seed

    void         seedChi2            (DeltaSeed* Seed, float Xc, float Yc, float& Chi2Par, float& Chi2Perp);
    //    void         seedDeltaChi2       (DeltaSeed* Seed, DeltaCandidate* Delta, float& Chi2Par, float& Chi2Perp);
  };
}
#endif
