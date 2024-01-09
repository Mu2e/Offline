#ifndef GeomPrimitives_Hole_HH
#define GeomPrimitives_Hole_HH
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
    Hole(){} // To build a Hole with the needed info.
    double               radius_;
    double               halfLength_;
    CLHEP::Hep3Vector         centerInMother_;
    std::string                orientation_;

  }; // end class Hole definition

} // end namespace mu2e

#endif //  GeomPrimitives_Hole_HH
