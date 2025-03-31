#ifndef RecoDataProducts_KKSHFlag_hh
#define RecoDataProducts_KKSHFlag_hh
//
// Bitflag used in KKStrawHit
//
//
// Original author David Brown
//
// Mu2e includes
#include "Offline/GeneralUtilities/inc/BitMap.hh"
#include <string>
#include <map>
namespace mu2e {

  struct KKSHFlagDetail {
    typedef unsigned mask_type;
    enum bit_type { tot=0, absdrift=1, driftdt=2,// if set, these values are used to constrain t0
      nhdrift=8, // if set, use drift radius to set null hit variance (otherwise use straw radius)
      longval=9, // if set, use time division to constrain longitudinal position
      annprob=15, // if set, interpret sign ANN probabilistically
      added=20 // record if the hit was added (otherwise original)
    };
    // functions needed for the BitMap template
    static std::string const& typeName();
    static std::map<std::string,mask_type> const& bitNames();
    static mask_type bit_to_mask( bit_type b){ return 1<<b; }
  };
  typedef BitMap<KKSHFlagDetail> KKSHFlag;
  typedef std::vector<mu2e::KKSHFlag> KKSHFlagCollection;
}
#endif
