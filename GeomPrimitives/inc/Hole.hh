#ifndef GeomPrimitives_Hole_HH
#define GeomPrimitives_Hole_HH
// This is the header file for the Notch class.  This class is a representation
// of a notch that might appear in any other geometry object.  It is
// simply a vector with the half sizes of a box, a Hep3Vector with the
// location of the center of the notch relative to the mother volume
// and the orientation of the notch as a string using the convention
// of the OrientationResolver.

//taken from Notch.hh for Holes

#include "CLHEP/Vector/ThreeVector.h"
#include <string>
#include <vector>

namespace mu2e {

  class Hole{
  public:
    Hole(double& rad, double& halfLen, CLHEP::Hep3Vector& location, std::string& ori)
      : radius_(rad),
        halfLength_ ( halfLen ),
        centerInMother_ ( location ),
        orientation_ ( ori )
    {}

    const double& getRad()       const { return radius_; }
    const double& getHalfLen()       const { return halfLength_; }
    const CLHEP::Hep3Vector&   getCenter()     const { return centerInMother_;}
    const std::string&         getOrient()     const { return orientation_; }

  private:
    Hole(){} // Have to build a Hole with the needed info.
    double               radius_;
    double               halfLength_;
    CLHEP::Hep3Vector         centerInMother_;
    std::string                orientation_;

  }; // end class Hole definition

} // end namespace mu2e

#endif //  GeomPrimitives_Hole_HH
