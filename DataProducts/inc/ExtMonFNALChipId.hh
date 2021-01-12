#ifndef DataProducts_ExtMonFNALChipId_hh
#define DataProducts_ExtMonFNALChipId_hh

// Identifier of a silicon chip in Mu2e ExtMonFNAL detector.
//
//
// Original author Andrei Gaponenko

#include <ostream>

#include "DataProducts/inc/ExtMonFNALModuleId.hh"

namespace mu2e {

  class ExtMonFNALChipId {
  public:

    ExtMonFNALChipId(const ExtMonFNALModuleId& module, unsigned int chipCol, unsigned int chipRow);

    // Default constructor should not be used by Mu2e code, but it is required by ROOT persistency
    ExtMonFNALChipId() : module_(), chipCol_(), chipRow_() {}

    const ExtMonFNALModuleId& module() const { return module_; }
    unsigned int chipCol() const { return chipCol_; }
    unsigned int chipRow() const { return chipRow_; }

    bool operator==( ExtMonFNALChipId const& rhs) const{
      return (module_ == rhs.module_)&&(chipRow_ == rhs.chipRow_)&&(chipCol_ == rhs.chipCol_);
    }

    bool operator!=( ExtMonFNALChipId const& rhs) const{
      return !(*this == rhs);
    }

    bool operator<( ExtMonFNALChipId const& rhs) const{
      return
        (module_ < rhs.module_) ||
        ((module_ == rhs.module_) && ((chipRow_ < rhs.chipRow_) ||
                                      ((chipRow_ == rhs.chipRow_) && (chipCol_ < rhs.chipCol_))
                                      )
         );
    }

  private:
    ExtMonFNALModuleId module_;
    unsigned int chipCol_;
    unsigned int chipRow_;
  };

  std::ostream& operator<<( std::ostream& os, const ExtMonFNALChipId& id);

}
#endif /* DataProducts_ExtMonFNALChipId_hh */
