//
//  Define the foils and bounding surfaces of the target
//  original author: David Brown (LBN) 2023
//
#ifndef KKGeom_StoppingTarget_hh
#define KKGeom_StoppingTarget_hh
#include "KinKal/Geometry/Cylinder.hh"
#include "KinKal/Geometry/Annulus.hh"
#include <vector>
namespace mu2e {
  namespace KKGeom {
    class StoppingTarget {
      public:
        // default constructor with nominal geometry
        StoppingTarget();
        // accessors
        auto const& outerCylinder() const { return outercyl_; }
        auto const& innerCylinder() const { return innercyl_; }
        auto const& foils() const { return foils_; }

      private:
        KinKal::Cylinder outercyl_; // outer cylinder
        KinKal::Cylinder innercyl_; //  inner cylinder
        std::vector<KinKal::Annulus> foils_; // target foils
    };
  }
}

#endif
