#ifndef RecoDataProducts_StrawHitFlag_hh
#define RecoDataProducts_StrawHitFlag_hh
//
// Class to describe flag bits used for straw hits
//
//
// Original author David Brown
//
// Mu2e includes
#include "Offline/GeneralUtilities/inc/BitMap.hh"
#include <string>
#include <map>
namespace mu2e {

  struct StrawHitFlagDetail {
// I need 32 bits for this class
    typedef unsigned mask_type;
    enum bit_type {stereo=0, energysel=1, radsel=2, timesel=3, bkgclust=5, bkg=6, isolated=7, outlier=8, other=9,
    tdiv=10, tclust=11,
    calosel=12, strawxtalk=13, elecxtalk=14, trksel=15,
    active=16,doca=17, resolvedphi=18,
    calopresel=19, intime=20, panelcombo=21, track=22,
    dead=23, noisy=24, nhitsel=25, sline=26};
    // functions needed for the BitMap template
    static std::string const& typeName();
    static std::map<std::string,mask_type> const& bitNames();
// maximum track Id I can flag
    static unsigned _maxTrkId;
    static mask_type bit_to_mask( bit_type b){ return 1<<b; }
  };
  typedef BitMap<StrawHitFlagDetail> StrawHitFlag;
  typedef std::vector<mu2e::StrawHitFlag> StrawHitFlagCollection;
}
#endif

