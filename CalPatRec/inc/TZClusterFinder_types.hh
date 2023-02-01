#ifndef __TrkPatRec_TZClusterFinder_types_hh__
#define __TrkPatRec_TZClusterFinder_types_hh__

#include <vector>

namespace art {
  class Event;
};

namespace fhicl {
  class ParameterSet;
};

class TH1D;

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Sequence.h"

#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitPosition.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/IntensityInfoTimeCluster.hh"
#include "Offline/Mu2eUtilities/inc/LsqSums4.hh"
#include "Offline/Mu2eUtilities/inc/McUtilsToolBase.hh"

namespace mu2e {

  namespace TZClusterFinderTypes {

//-----------------------------------------------------------------------------
// data structure shared by TZClusterFinder with its plugins
//-----------------------------------------------------------------------------

    struct Config {
      fhicl::Atom<std::string> tool_type             {fhicl::Name("tool_type"             ), fhicl::Comment("tool type: TZClusterFinderDiag") };
      fhicl::Atom<int>         mcTruth               {fhicl::Name("mcTruth"               ), fhicl::Comment("MC truth")                       };
      fhicl::Atom<int>         simIDThresh           {fhicl::Name("simIDThresh"           ), fhicl::Comment("number of straw hits")           };

      fhicl::Table<McUtilsToolBase::Config> mcUtils{fhicl::Name("mcUtils"       ), fhicl::Comment("MC Diag plugin") };
    };


    struct cHit {
      cHit():hIndex(-1),hTime(0.),hZpos(0.),nStrawHits(0),hIsUsed(0){}
      int    hIndex;
      double hTime;
      double hZpos;
      int    nStrawHits;
      int    hIsUsed;
    };

    struct plnData {
      std::vector<cHit> plnHits;
    };

    struct chunkInfo {
      std::vector<int> hIndices;
      double avgTime;
      double avgZpos;
      int nHits; // combo hits
      int nStrawHits;
      int nrgSelection; // 1 if passes energy selection (CE), 0 if not (protons)
      int nCombines;
      int caloIndex = -1;
      double lSlope;
      double lIntercept;
      double chi2Dof;
    };

    struct Data_t {

      const art::Event*               _event;
      const TimeCluster*              _timeCluster;
      const ComboHitCollection*       chcol;
      const ComboHitCollection*       chcol2; // for tool to get simID info before selection cuts
      const CaloClusterCollection*    ccCollection;
      const StrawHitFlagCollection*   shfcol;
      TimeClusterCollection*          _tcColl;       // 'tcColl': time cluster collection
      IntensityInfoTimeCluster*       _iiTC;
      int                             _nTZClusters;

      // members that may eventually be moved to TimeCluster.hh
      std::vector<double> lineSlope;
      std::vector<double> lineIntercept;
      std::vector<double> chi2DOF;

      // clear function for members that need to be cleared before processing new event
      void clearRelevant() {
        lineSlope.clear();
        lineIntercept.clear();
        chi2DOF.clear();
      }
    };

    struct pair {
      int first;
      int second;
    };

    // create data structure for variables and objects that will facilitate cluster finding
    struct facilitateVars {

      std::array<plnData, StrawId::_nplanes> cHits;
      std::vector<chunkInfo> chunks;
      chunkInfo _chunkInfo;
      TimeCluster _clusterInfo;
      ::LsqSums4 fitHits;
      std::vector<pair> holdIndices;
      int seedIndice;
      double seedTime;
      double seedZpos;
      int seedNRGselection;
      pair indicePair;
      int startIndex;
      int testIndice;
      double testTime;
      double testZpos;
      int testNRGselection;
      int nHitsInChunk;
      int nStrawHitsInChunk;
      double totalTime;
      double totalZpos;
      bool moreCombines;
      int chunkOneIdx;
      int chunkTwoIdx;
      double chi2seed;
      double chi2combineTest;
      double biggestChi2combine;

      void clear_cHits() { for (size_t i=0; i<cHits.size(); i++) cHits[i].plnHits.clear(); }
      void clear_clusterInfo() { _clusterInfo._strawHitIdxs.clear(); }
      void clear_chunkInfo() { _chunkInfo.hIndices.clear(); }
      void clear_chunks() { chunks.clear(); }

    };

    struct mcSimIDs {
      int simID;
      int pdgID;
      int nHits; // number of straw hits corresponding to simID
    };
  }
}
#endif
