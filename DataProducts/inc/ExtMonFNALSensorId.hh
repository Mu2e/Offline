#ifndef TrackerGeom_ExtMonFNALSensorId_hh
#define TrackerGeom_ExtMonFNALSensorId_hh

// Identifier of a silicon sensor wafer.  At the moment it is just an
// integer (the plane number), but this may change if we decide to
// have more than one sensor per plane.
//
// $Id: ExtMonFNALSensorId.hh,v 1.3 2012/08/28 05:02:51 gandr Exp $
// $Author: gandr $
// $Date: 2012/08/28 05:02:51 $
//
// Original author Andrei Gaponenko

#include <ostream>

namespace mu2e {

  class ExtMonFNALSensorId{
  public:

    static const unsigned int NOPLANE = -1u;

    // No automatic conversion of int to ExtMonFNALSensorId.
    explicit ExtMonFNALSensorId(unsigned int plane) : plane_(plane) {}

    // zero based
    unsigned int plane() const { return plane_;}

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
    unsigned int plane_;
  };

  inline std::ostream& operator<<( std::ostream& os, const ExtMonFNALSensorId& id) {
    return os<<id.plane();
  }
}
#endif /* TrackerGeom_ExtMonFNALSensorId_hh */
