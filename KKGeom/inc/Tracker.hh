//
//  Define the nominal tracker boundary and reference surfaces, used to extrapolate and sample KinKal track fits, and to build
//  the passive materials in the fit
//  original author: David Brown (LBN) 2023
//
#ifndef KKGeom_Tracker_hh
#define KKGeom_Tracker_hh
#include "KinKal/Geometry/Cylinder.hh"
#include "KinKal/Geometry/Disk.hh"
#include "KinKal/Geometry/Annulus.hh"
#include <exception>
namespace mu2e {
  namespace KKGeom {
    class Tracker {
      public:
        // default constructor with nominal geometry
        Tracker();
        // accessors
        // define the positions corresponding to the G4 virtual detectors
        // 'entrance' refers to upstream (negative z) end, 'middle' to the middle of the tracker, 'exit' refers to downstream (positive z) end
        auto const& outerCylinder() const { return outercyl_; }
        auto const& innerCylinder() const { return innercyl_; }
        auto const& entDisk() const { return ent_; }
        auto const& midDisk() const { return mid_; }
        auto const& exitDisk() const { return exit_; }

      private:
        KinKal::Cylinder outercyl_; // outer cylinder
        KinKal::Cylinder innercyl_; //  inner cylinder
        KinKal::Disk ent_, mid_, exit_; // standard reference point planes
        // TODO: add passive material for manifolds, absorber rings, beams, ...
 //       std::vector<Annulus> manifolds_; // manifolds
    };
  }
}

#endif
