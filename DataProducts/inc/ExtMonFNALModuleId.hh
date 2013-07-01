#ifndef TrackerGeom_ExtMonFNALModuleId_hh
#define TrackerGeom_ExtMonFNALModuleId_hh

// Modified ExtMonFNALSensorId.hh (original author Andrei Gaponenko)
// Author Evan Schiewe

#include <ostream>

namespace mu2e {

  class ExtMonFNALModuleId{
  public:
    enum Side {FRONT, BACK};  // FRONT is upstream, BACK is downstream
    enum Rotation {UP, DOWN}; // direction refers to placement of readout cables (ie. UP means end of module with cables points up)
    static const unsigned int NOPLANE = -1u;

    explicit ExtMonFNALModuleId(unsigned int plane, Side side, Rotation rot, unsigned int module) : plane_(plane), side_(side), rot_(rot), module_(module) {}

    // zero based
    unsigned int plane() const { return plane_;}
    Side side() const { return side_; }
    Rotation rot() const { return rot_; }
    unsigned int module() const { return module_;}

    bool operator==( ExtMonFNALModuleId const& rhs) const{
      return (plane_ == rhs.plane_ && side_ == rhs.side_ && rot_ == rhs.rot_ && module_ == rhs.module_);
    }

    bool operator!=( ExtMonFNALModuleId const& rhs) const{
      return !(*this == rhs);
    }

    bool operator<( ExtMonFNALModuleId const& rhs) const{
      return ((plane_ < rhs.plane_) || 
              ((plane_ == rhs.plane_) && ((side_ == BACK) && (rhs.side_ == FRONT))) || 
              ((plane_ == rhs.plane_) && (side_ == rhs.side_) && (module_ < rhs.module_)));
    }

    // Default constructor is required by ROOT persistency
    ExtMonFNALModuleId() : plane_(NOPLANE), side_(FRONT), rot_(UP), module_() {}

  private:
    unsigned int plane_;
    Side side_;
    Rotation rot_;
    unsigned int module_;
  };

  inline std::ostream& operator<<( std::ostream& os, const ExtMonFNALModuleId& id) {
    enum Side {FRONT, BACK};
    std::string side;
    switch(id.side()) {
    case (FRONT): side = "FRONT"; break;
    case (BACK): side = "BACK"; break;
    default: side = "UNKNOWN"; break;
    }

    enum Rotation {UP, DOWN};
    std::string rot;
    switch(id.rot()) {
    case (UP): rot = "UP"; break;
    case (DOWN): rot = "DOWN"; break;
    default: rot = "UNKNOWN"; break;
    }
    return os<<id.plane()
             <<","<<side
             <<","<<rot
             <<","<<id.module();
  }
}
#endif /* TrackerGeom_ExtMonFNALModuleId_hh */
