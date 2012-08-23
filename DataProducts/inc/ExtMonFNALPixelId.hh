#ifndef TrackerGeom_ExtMonFNALPixelId_hh
#define TrackerGeom_ExtMonFNALPixelId_hh

// Identifier of a silicon pixel in Mu2e ExtMonFNAL detector.
//
// $Id: ExtMonFNALPixelId.hh,v 1.1 2012/08/23 23:41:34 gandr Exp $
// $Author: gandr $
// $Date: 2012/08/23 23:41:34 $
//
// Original author Andrei Gaponenko

#include <ostream>

#include "DataProducts/inc/ExtMonFNALChipId.hh"

namespace mu2e {

  class ExtMonFNALPixelId {
  public:

    ExtMonFNALPixelId(const ExtMonFNALChipId& chip, int col, int row);

    // Default constructor should not be used by Mu2e code, but it is required by ROOT persistency
    ExtMonFNALPixelId() : chip_(), col_(), row_() {}

    const ExtMonFNALChipId& chip() const { return chip_; }
    int   col() const { return col_; }
    int   row() const { return row_; }

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
    int col_;
    int row_;
  };

  std::ostream& operator<<( std::ostream& os, const ExtMonFNALPixelId& id);

}
#endif /* TrackerGeom_ExtMonFNALPixelId_hh */
