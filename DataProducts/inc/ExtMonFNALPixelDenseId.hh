#ifndef TrackerGeom_ExtMonFNALPixelDenseId_hh
#define TrackerGeom_ExtMonFNALPixelDenseId_hh

// Sequential number of a silicon pixel in Mu2e ExtMonFNAL detector.
// Zero based.
//
// $Id: ExtMonFNALPixelDenseId.hh,v 1.1 2012/08/28 05:05:27 gandr Exp $
// $Author: gandr $
// $Date: 2012/08/28 05:05:27 $
//
// Original author Andrei Gaponenko

#include <ostream>

namespace mu2e {

  class ExtMonFNALPixelDenseId {
  public:

    static const unsigned int NOPIXEL = -1u;

    explicit ExtMonFNALPixelDenseId(unsigned int pix = NOPIXEL) : pix_(pix) {}

    bool operator==( ExtMonFNALPixelDenseId const& rhs) const{
      return pix_ == rhs.pix_;
    }

    bool operator!=( ExtMonFNALPixelDenseId const& rhs) const{
      return !(*this == rhs);
    }

    bool operator<( ExtMonFNALPixelDenseId const& rhs) const{
      return pix_ < rhs.pix_;
    }

    unsigned int number() const { return pix_; }

  private:
    unsigned int pix_;
  };

  inline std::ostream& operator<<( std::ostream& os, const ExtMonFNALPixelDenseId& id) {
    return os<<id.number();
  }

}
#endif /* TrackerGeom_ExtMonFNALPixelDenseId_hh */
