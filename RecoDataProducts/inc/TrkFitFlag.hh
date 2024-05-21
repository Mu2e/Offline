#ifndef RecoDataProducts_TrkFitFlag_hh
#define RecoDataProducts_TrkFitFlag_hh
//
// Class to describe flag bits used for track fits
//
//
// Original author David Brown
//
// Mu2e includes
#include "Offline/GeneralUtilities/inc/BitMap.hh"
#include <string>
#include <map>
namespace mu2e {

  struct TrkFitFlagDetail {
    // I need 32 bits for this class
    typedef unsigned mask_type;
    enum bit_type {hitsOK=0,circleOK,phizOK,helixOK,seedOK,kalmanOK,circleInit,phizInit,
    circleConverged,phizConverged,helixConverged,seedConverged,kalmanConverged,
    MatCorr, BFCorr, FitOK,
    KSF=16, KFF, TPRHelix, CPRHelix, Straight, KKLoopHelix,KKCentralHelix,KKLine, APRHelix, MPRHelix, MCSeed=31};

    // functions needed for the BitMap template
    static std::string const& typeName();
    static std::map<std::string,mask_type> const& bitNames();
    static mask_type bit_to_mask( bit_type b){ return 1<<b; }
  };
  typedef BitMap<TrkFitFlagDetail> TrkFitFlag;

}
#endif

