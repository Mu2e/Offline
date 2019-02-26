
#ifndef DataProducts_ExtMonFNALPlaneId_hh
#define DataProducts_ExtMonFNALPlaneId_hh

#include <ostream>

// currently just an integer for each of the 8 planes.

namespace mu2e {

  class ExtMonFNALPlaneId{
  public:

    static const unsigned int NOPLANE = -1u;

    // No automatic conversion of int to ExtMonFNALPlaneId.
    explicit ExtMonFNALPlaneId(unsigned int plane) : plane_(plane) {}

    // zero based
    unsigned int plane() const { return plane_;}

    bool operator==( ExtMonFNALPlaneId const& rhs) const{
      return (plane_ == rhs.plane_);
    }

    bool operator!=( ExtMonFNALPlaneId const& rhs) const{
      return !(*this == rhs);
    }

    bool operator<( ExtMonFNALPlaneId const& rhs) const{
      return (plane_ < rhs.plane_);
    }

    // Default constructor is required by ROOT persistency
    ExtMonFNALPlaneId() : plane_(NOPLANE) {}

  private:
    unsigned int plane_;
  };

  inline std::ostream& operator<<( std::ostream& os, const ExtMonFNALPlaneId& id) {
    return os<<id.plane();
  }
}
#endif /* DataProducts_ExtMonFNALPlaneId_hh */
