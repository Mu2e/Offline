#ifndef DataProducts_ExtMonFNALModuleId_hh
#define DataProducts_ExtMonFNALModuleId_hh

// Author Evan Schiewe

#include <ostream>

namespace mu2e {

  class ExtMonFNALModuleId{
  public:
    static const unsigned int NOPLANE = -1u;

    explicit ExtMonFNALModuleId(unsigned int plane,  unsigned int module) : plane_(plane), module_(module) {}
      
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
#endif /* DataProducts_ExtMonFNALModuleId_hh */
