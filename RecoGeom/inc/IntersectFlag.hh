#ifndef RecoGeom_IntersectFlag_hh
#define RecoGeom_IntersectFlag_hh
//
// describe flag bits used for reco intersection
//
//
// Original author David Brown
//
// Mu2e includes
#include "Offline/GeneralUtilities/inc/BitMap.hh"
#include <string>
#include <map>
namespace mu2e {
  namespace RecoGeom {
    struct IntersectFlagDetail {
      // I need 16 bits for this class
      typedef unsigned short mask_type;
      enum bit_type {onsurface=0,inbounds=1 };
      // functions needed for the BitMap template
      static std::string const& typeName();
      static std::map<std::string,mask_type> const& bitNames();
      // maximum track Id I can flag
      static unsigned _maxTrkId;
      static mask_type bit_to_mask( bit_type b){ return 1<<b; }
    };
    typedef BitMap<IntersectFlagDetail> IntersectFlag;
  }
}
#endif
