#ifndef CalPatRec_AgnosticHelixFinder_types_hh
#define CalPatRec_AgnosticHelixFinder_types_hh

#include "fhiclcpp/types/Atom.h"
#include "art/Framework/Principal/Event.h"

#include "Offline/Mu2eUtilities/inc/LsqSums2.hh"
#include "Offline/Mu2eUtilities/inc/LsqSums4.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"

namespace mu2e {

  namespace AgnosticHelixFinderTypes {

    //-----------------------------------------------------------------------------
    // Used in the AgnosticHelixFinder
    //-----------------------------------------------------------------------------

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

    struct lineInfo {
      float                 zMin;
      float                 zMax;
      ::LsqSums2            fitter;
      std::vector<size_t>   tcHitsIndices;
      std::vector<int>      helixPhiCorrections; // integer for 2pi ambiguity
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

    //-----------------------------------------------------------------------------
    // Used in the AgnosticHelixFinderDiag
    //-----------------------------------------------------------------------------

    struct Config {
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<std::string> tool_type{Name("tool_type")      , Comment("tool type: AgnosticHelixFinderDiag")};
      fhicl::Atom<std::string> simTag   {Name("simTag")         , Comment("Sim particle collection tag"), ""};
      fhicl::Atom<std::string> digTag   {Name("strawDigiMCTag") , Comment("Straw digi MC collection tag"), ""};
      fhicl::Atom<bool>        display  {Name("display")        , Comment("Run the display"), false};
      fhicl::Atom<int>         debug    {Name("debug")          , Comment("Tool debug printouts"), 0};
    };

    struct tripletInfo { // for diagnostics
      const triplet* trip;
      float   radius = 0.f;
      float   xC = 0.f;
      float   yC = 0.f;
    };

    struct tcInfo {
      int     nHelices;
      int     nComboHits;
      int     nStrawHits;
    };

    struct hsInfo {
      float eDepAvg;
      float bz0;
      const HelixSeed* seed;
      hsInfo(float EDep = 0., float B = 0., const HelixSeed* Seed = nullptr) :
        eDepAvg(EDep), bz0(B), seed(Seed) {}
      hsInfo(float B, const HelixSeed* Seed) : hsInfo(0., B, Seed) {}
    };

    struct lineSegmentInfo {
      float   chi2dof;
      float   maxHitGap;
    };

    struct finalLineInfo{
      float nHitsRatio;
    };

    struct diagInfo {
      int                            nHelices;
      int                            nTimeClusters;
      std::vector<tcInfo>            timeClusterData;
      std::vector<hsInfo>            helixSeedData;
      std::vector<lineSegmentInfo>   lineSegmentData;
      std::vector<finalLineInfo>     lineInfoData;

      // pointer to the art event
      const art::Event* event = nullptr;

      // pointers to data currently under evaluation
      const HelixSeed* hseed = nullptr;
      const std::vector<cHit>* tcHits = nullptr;
      const std::vector<lineInfo>* seedPhiLines = nullptr;
      const ComboHitCollection* chColl = nullptr;

      XYZVectorF caloPos;
      XYZVectorF targPos;
      tripletInfo tripInfo;
      float bz0 = 0.f;

      int diagLevel = 0;
    };

    // Special "hits" used in the helix finding
    enum HitType {
      CALOCLUSTER = -1,
      STOPPINGTARGET = -2
    };

    // Enums for the diagnostic tool
    enum DIAG {kBegin = 1,
      kTriplet, kCircle, kSegments, kLine, kRecover, kViability, kHelix, kFinal, // helix fit stages
      kTimeCluster, kEnd};

  } // namespace AgnosticHelixFinderTypes
} // namespace mu2e
#endif
