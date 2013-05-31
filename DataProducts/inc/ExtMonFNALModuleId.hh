#ifndef TrackerGeom_ExtMonFNALModuleId_hh
#define TrackerGeom_ExtMonFNALModuleId_hh

// Modified ExtMonFNALSensorId.hh (original author Andrei Gaponenko)
// Author Evan Schiewe

#include <ostream>

namespace mu2e {

  class ExtMonFNALModuleId{
  public:
    enum Side {FRONT, BACK}; // not using "enum class" because of scoping issues with L46 (inline ostream)
    static const unsigned int NOPLANE = -1u;

    explicit ExtMonFNALModuleId(unsigned int plane, Side side, unsigned int module) : plane_(plane), side_(side), module_(module) {}

    // zero based
    unsigned int plane() const { return plane_;}
    Side side() const { return side_; }
    unsigned int module() const { return module_;}

    bool operator==( ExtMonFNALModuleId const& rhs) const{
      return (plane_ == rhs.plane_ && side_ == rhs.side_ && module_ == rhs.module_);
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
    ExtMonFNALModuleId() : plane_(NOPLANE), side_(FRONT), module_() {}

  private:
    unsigned int plane_;
    Side side_;
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
    return os<<id.plane()<<","
             <<side
             <<","<<id.module();
  }
}
#endif /* TrackerGeom_ExtMonFNALModuleId_hh */
