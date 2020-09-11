#ifndef DataProducts_ExtMonFNALPixelId_hh
#define DataProducts_ExtMonFNALPixelId_hh

// Identifier of a silicon pixel in Mu2e ExtMonFNAL detector.
//
//
// Original author Andrei Gaponenko

#include <ostream>

#include "DataProducts/inc/ExtMonFNALChipId.hh"

namespace mu2e {

  class ExtMonFNALPixelId {
  public:

    ExtMonFNALPixelId(const ExtMonFNALChipId& chip, unsigned int col, unsigned int row);

    // Default constructor should not be used by Mu2e code, but it is required by ROOT persistency
    ExtMonFNALPixelId() : chip_(), col_(), row_() {}

    const ExtMonFNALChipId& chip() const { return chip_; }

    // Row and column numbers are zero based - this is not the hardware numbering
    unsigned int col() const { return col_; }
    unsigned int row() const { return row_; }

    bool operator==( ExtMonFNALPixelId const& rhs) const{
      return (chip_ == rhs.chip_)&&(row_ == rhs.row_)&&(col_ == rhs.col_);
    }

    bool operator!=( ExtMonFNALPixelId const& rhs) const{
      return !(*this == rhs);
    }

    bool operator<( ExtMonFNALPixelId const& rhs) const{
      return
        (chip_ < rhs.chip_) ||
        ((chip_ == rhs.chip_) && ((row_ < rhs.row_) ||
                                  ((row_ == rhs.row_) && (col_ < rhs.col_))
                                  )
         );
    }

  private:
    ExtMonFNALChipId chip_;
    unsigned int col_;
    unsigned int row_;
  };

  std::ostream& operator<<( std::ostream& os, const ExtMonFNALPixelId& id);

}
#endif /* DataProducts_ExtMonFNALPixelId_hh */
