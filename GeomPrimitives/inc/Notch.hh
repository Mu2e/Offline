#ifndef GeomPrimitives_Notch_HH
#define GeomPrimitives_Notch_HH
// This is the header file for the Notch class.  This class is a representation
// of a notch that might appear in any other geometry object.  It is
// simply a vector with the half sizes of a box, a Hep3Vector with the 
// location of the center of the notch relative to the mother volume
// and the orientation of the notch as a string using the convention
// of the OrientationResolver.

#include "CLHEP/Vector/ThreeVector.h"
#include <string>
#include <vector>

namespace mu2e {

  class Notch{
  public:
    Notch( std::vector<double>& halfSizes, CLHEP::Hep3Vector& location,
	   std::string& orientString )
      : halfDims_ ( halfSizes ),
	centerInMother_ ( location ),
	orientation_ ( orientString )
    {}

    const std::vector<double>& getDims()       const { return halfDims_; }
    const CLHEP::Hep3Vector&   getCenter()     const { return centerInMother_;}
    const std::string&         getOrient()     const { return orientation_; }

  private:
    Notch(){} // Have to build a Notch with the needed info.
    std::vector<double>           halfDims_;
    CLHEP::Hep3Vector             centerInMother_;
    std::string                   orientation_;

  }; // end class Notch definition

} // end namespace mu2e

#endif //  GeomPrimitives_Notch_HH

