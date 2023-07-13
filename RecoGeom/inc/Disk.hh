//
//  Description of an annular planar section
//  original author: David Brown (LBN) 2023
//
#ifndef RecoGeom_Disk_hh
#define RecoGeom_Disk_hh
#include "Offline/RecoGeom/inc/Annulus.hh"
namespace mu2e {
  namespace RecoGeom {
    class Disk : public Annulus {
      public:
        // construct from necessary parameters
        Disk(XYZVectorD const& norm, XYZVectorD const& center, double radius) : Annulus(norm,center,0.0,radius) {}
    };
  }
}
std::ostream& operator <<(std::ostream& ost, mu2e::RecoGeom::Disk const& disk);
#endif
