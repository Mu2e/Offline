#ifndef TrackerGeom_ExtMonFNALModuleId_hh
#define TrackerGeom_ExtMonFNALModuleId_hh

// Author Evan Schiewe

#include <ostream>

namespace mu2e {

  class ExtMonFNALModuleId{
  public:
    static const unsigned int NOPLANE = -1u;

    explicit ExtMonFNALModuleId(unsigned int plane,  unsigned int module) : plane_(plane), module_(module) {}
    explicit ExtMonFNALModuleId(unsigned int mod)
    {
      //TODO: remove this constructor and fix ExtMonFNALPixelIdConverter
      plane_ = (mod < 12 ? (mod/3) : (mod/2 - 2));
      module_ = (plane_ < 4 ? (mod-plane_*3) : (mod -12 - ((plane_-4)*2)));
    }
      
    // zero based
    unsigned int plane() const { return plane_;}
    unsigned int number() const { return module_;}

    bool operator==( ExtMonFNALModuleId const& rhs) const{
      return (plane_ == rhs.plane_ && module_ == rhs.module_);
    }

    bool operator!=( ExtMonFNALModuleId const& rhs) const{
      return !(*this == rhs);
    }

    bool operator<( ExtMonFNALModuleId const& rhs) const{
      return ((plane_ < rhs.plane_) || 
              ((plane_ == rhs.plane_) && (module_ < rhs.module_)));
    }

    // Default constructor is required by ROOT persistency
    ExtMonFNALModuleId() : plane_(NOPLANE),  module_() {}
    
  private:
    unsigned int plane_;
    unsigned int module_;
  };
  
  inline std::ostream& operator<<( std::ostream& os, const ExtMonFNALModuleId& id) {
    return os<<id.plane()
             <<","<<id.number();
  }
}
#endif /* TrackerGeom_ExtMonFNALModuleId_hh */
