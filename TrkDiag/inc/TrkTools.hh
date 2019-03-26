//
// Namespace for collecting tools used in TrkDiag tree filling
// Original author: A. Edmonds (November 2018)
//
#ifndef TrkDiag_TrkTools_hh
#define TrkDiag_TrkTools_hh
#include "RecoDataProducts/inc/StrawHitIndex.hh"
#include "RecoDataProducts/inc/KalSeed.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "RecoDataProducts/inc/TrkQual.hh"
#include "RecoDataProducts/inc/RecoCount.hh"
//
// Namespace for collecting tools used in TrkDiag tree filling
// Original author: A. Edmonds (November 2018)
//
#ifndef TrkDiag_TrkTools_hh
#define TrkDiag_TrkTools_hh
#include "RecoDataProducts/inc/StrawHitIndex.hh"
#include "RecoDataProducts/inc/KalSeed.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "RecoDataProducts/inc/TrkQual.hh"
#include "RecoDataProducts/inc/RecoCount.hh"

#include "CalorimeterGeom/inc/DiskCalorimeter.hh"

#include "TrkDiag/inc/HitCount.hh"
#include "TrkDiag/inc/TrkInfo.hh"
#include "TrkDiag/inc/TrkStrawHitInfo.hh"
#include "TrkDiag/inc/TrkStrawMatInfo.hh"
#include "TrkDiag/inc/TrkCaloHitInfo.hh"
#include "TrkDiag/inc/TrkQualInfo.hh"
#include "TrkDiag/inc/HelixInfo.hh"

#include <vector>
#include <functional>
namespace mu2e {
  namespace TrkTools {
    // count different hit types
    void countHits(const std::vector<TrkStrawHitSeed>& hits, unsigned& nhits, unsigned& nactive, unsigned& ndouble, unsigned& ndactive, unsigned& nnullambig);
    void countStraws(const std::vector<TrkStraw>& straws, unsigned& nmat, unsigned& nmatactive, double& radlen);

    // fill various Info structs
    void fillHitCount(const StrawHitFlagCollection& flags, HitCount& hitcount);
    void fillHitCount(RecoCount const& nrec, HitCount& hitcount);
    void fillTrkInfo(const KalSeed& kseed,TrkInfo& trkinfo);

    void fillHitInfo(const KalSeed& kseed, std::vector<TrkStrawHitInfo>& tshinfos );
    void fillMatInfo(const KalSeed& kseed, std::vector<TrkStrawMatInfo>& tminfos );
    void fillCaloHitInfo(const KalSeed& kseed, Calorimeter const& calo, TrkCaloHitInfo& tchinfo );
    void fillTrkQualInfo(const TrkQual& tqual, TrkQualInfo& trkqualInfo);
    void fillHelixInfo(const KalSeed& kseed, double bz0, HelixInfo& hinfo); 
  }
}

#endif

#include "CalorimeterGeom/inc/DiskCalorimeter.hh"

#include "TrkDiag/inc/HitCount.hh"
#include "TrkDiag/inc/TrkInfo.hh"
#include "TrkDiag/inc/TrkStrawHitInfo.hh"
#include "TrkDiag/inc/TrkStrawMatInfo.hh"
#include "TrkDiag/inc/TrkCaloHitInfo.hh"
#include "TrkDiag/inc/TrkQualInfo.hh"
#include "TrkDiag/inc/HelixInfo.hh"

#include <vector>
#include <functional>
namespace mu2e {
  namespace TrkTools {
    // count different hit types
    void countHits(const std::vector<TrkStrawHitSeed>& hits, unsigned& nhits, unsigned& nactive, unsigned& ndouble, unsigned& ndactive, unsigned& nnullambig);
    void countStraws(const std::vector<TrkStraw>& straws, unsigned& nmat, unsigned& nmatactive, double& radlen);

    // fill various Info structs
    void fillHitCount(const StrawHitFlagCollection& flags, HitCount& hitcount);
    void fillHitCount(RecoCount const& nrec, HitCount& hitcount);
    void fillTrkInfo(const KalSeed& kseed,TrkInfo& trkinfo);

    void fillHitInfo(const KalSeed& kseed, std::vector<TrkStrawHitInfo>& tshinfos );
    void fillMatInfo(const KalSeed& kseed, std::vector<TrkStrawMatInfo>& tminfos );
    void fillCaloHitInfo(const KalSeed& kseed, Calorimeter const& calo, TrkCaloHitInfo& tchinfo );
    void fillTrkQualInfo(const TrkQual& tqual, TrkQualInfo& trkqualInfo);
    void fillHelixInfo(const KalSeed& kseed, double bz0, HelixInfo& hinfo); 
  }
}

#endif
