#ifndef TrackerGeom_ExtMonFNALSensorId_hh
#define TrackerGeom_ExtMonFNALSensorId_hh

// Identifier of a silicon sensor wafer.  At the moment it is just an
// integer (the plane number), but this may change if we decide to
// have more than one sensor per plane.
//
// $Id: ExtMonFNALSensorId.hh,v 1.2 2012/08/23 23:41:34 gandr Exp $
// $Author: gandr $
// $Date: 2012/08/23 23:41:34 $
//
// Original author Andrei Gaponenko

#include <ostream>

namespace mu2e {

  class ExtMonFNALSensorId{
  public:

    static const int NOPLANE = -1;

    // No automatic conversion of int to ExtMonFNALSensorId.
    explicit ExtMonFNALSensorId(int plane) : plane_(plane) {}

    int plane() const { return plane_;}

    bool operator==( ExtMonFNALSensorId const& rhs) const{
      return (plane_ == rhs.plane_);
    }

    bool operator!=( ExtMonFNALSensorId const& rhs) const{
      return !(*this == rhs);
    }

    bool operator<( ExtMonFNALSensorId const& rhs) const{
      return (plane_ < rhs.plane_);
    }

    // Default constructor is required by ROOT persistency
    ExtMonFNALSensorId() : plane_(NOPLANE) {}

  private:
    int plane_;
  };

  inline std::ostream& operator<<( std::ostream& os, const ExtMonFNALSensorId& id) {
    return os<<id.plane();
  }
}
#endif /* TrackerGeom_ExtMonFNALSensorId_hh */
