#ifndef RecoGeom_IntersectFlag_hh
#define RecoGeom_IntersectFlag_hh
//
// describe flag bits used for reco intersection
//
//
// Original author David Brown
//
#include <string>
#include <map>
namespace mu2e {
  namespace RecoGeom {
    struct IntersectFlag {
      IntersectFlag() : onsurface_(false), inbounds_(false) {}
      bool onsurface_; // intersection is on the surface
      bool inbounds_;  // intersection is inside the surface boundaries
    };
  }
}
std::ostream& operator <<(std::ostream& ost, mu2e::RecoGeom::IntersectFlag const& iflag);
#endif
