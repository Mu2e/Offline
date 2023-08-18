#ifndef CalPatRec_TZClusterFinder_types_hh
#define CalPatRec_TZClusterFinder_types_hh

#include <vector>

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Sequence.h"

#include "art/Framework/Principal/Event.h"

#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitPosition.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/IntensityInfoTimeCluster.hh"
#include "Offline/Mu2eUtilities/inc/LsqSums2.hh"
#include "Offline/Mu2eUtilities/inc/McUtilsToolBase.hh"

namespace mu2e {

  namespace TZClusterFinderTypes {

//-----------------------------------------------------------------------------
// data structure shared by TZClusterFinder with its plugins
//-----------------------------------------------------------------------------

    struct Config {
      fhicl::Atom<std::string>              tool_type   {fhicl::Name("tool_type"  ), fhicl::Comment("tool type: TZClusterFinderDiag") };
      fhicl::Atom<int>                      mcTruth     {fhicl::Name("mcTruth"    ), fhicl::Comment("MC truth")                       };
      fhicl::Atom<int>                      simIDThresh {fhicl::Name("simIDThresh"), fhicl::Comment("number of straw hits")           };
      fhicl::Table<McUtilsToolBase::Config> mcUtils     {fhicl::Name("mcUtils"    ), fhicl::Comment("MC Diag plugin")                 };
    };


    struct cHit {
      int    hIndex;
      double hTime;
      double hWeight;
      double hZpos;
      int    nStrawHits;
      int    hIsUsed;
    };

    struct plnData {
      std::vector<cHit> plnHits;
    };

    struct chunkInfo {
      std::vector<int> hIndices;
      ::LsqSums2       fitter;
      double           avgTime;
      double           avgZpos;
      int              nHits; // combo hits
      int              nStrawHits;
      int              nrgSelection; // 1 if passes energy selection (CE), 0 if not (protons)
      int              nCombines;
      int              caloIndex;
    };

    struct Data_t {

      const art::Event*               _event;
      const TimeCluster*              _timeCluster;
      const ComboHitCollection*       _chColl;
      const ComboHitCollection*       _chColl2; // for tool to get simID info before selection cuts
      const CaloClusterCollection*    _ccColl;
      TimeClusterCollection*          _tcColl; // 'tcColl': time cluster collection
      IntensityInfoTimeCluster*       _iiTC;
      int                             _nTZClusters;

      // diagnostic data members used in TZ tool
      std::vector<double> lineSlope;
      std::vector<double> lineIntercept;
      std::vector<double> chi2DOF;

      // clear diagnostic information used in TZ tool
      void clearDiagInfo() {
        lineSlope.clear();
        lineIntercept.clear();
        chi2DOF.clear();
      }
    };

    struct cPair {
      int first;
      int second;
    };

    // create data structure for variables and objects that will facilitate cluster finding
    struct facilitateVars {

      std::array<plnData, StrawId::_nplanes> cHits;
      std::vector<chunkInfo>                 chunks;
      std::vector<cPair>                     holdIndices;
      chunkInfo   _chunkInfo;
      TimeCluster _clusterInfo;
      cPair       _indicePair;
      int         seedIndice;
      double      seedTime;
      double      seedWeight;
      double      seedZpos;
      int         seedNRGselection;
      int         startIndex;
      int         testIndice;
      double      testTime;
      double      testWeight;
      double      testZpos;
      int         testNRGselection;
      int         nHitsInChunk;
      int         nStrawHitsInChunk;
      double      totalTime;
      double      totalZpos;
      bool        moreCombines;

      void clear_cHits()       { for (size_t i=0; i<cHits.size(); i++) cHits[i].plnHits.clear(); }
      void clear_clusterInfo() { _clusterInfo._strawHitIdxs.clear();                             }
      void clear_chunkInfo()   { _chunkInfo.hIndices.clear(); _chunkInfo.fitter.clear();         }
      void clear_chunks()      { chunks.clear();                                                 }

    };

    struct mcSimIDs {
      int simID;
      int pdgID;
      int nHits; // number of straw hits corresponding to simID
    };
  }
}
#endif
