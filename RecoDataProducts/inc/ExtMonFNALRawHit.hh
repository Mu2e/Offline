#ifndef RecoDataProducts_ExtMonFNALRawHit_hh
#define RecoDataProducts_ExtMonFNALRawHit_hh
//
// Unpacked data from the ExtMonFNAL pixel detector.
//
// $Id: ExtMonFNALRawHit.hh,v 1.1 2012/08/23 23:41:35 gandr Exp $
// $Author: gandr $
// $Date: 2012/08/23 23:41:35 $
//
// Original author Andrei Gaponenko
//

#include <ostream>

#include "DataProducts/inc/ExtMonFNALPixelId.hh"

namespace mu2e {

  class ExtMonFNALRawHit {
  public:

    ExtMonFNALRawHit(const ExtMonFNALPixelId& pix,
                      int clock,
                      unsigned tot)
      : pixelId_(pix)
      , clock_(clock)
      , tot_(tot)
    {}

    // Default constructor should not be used by Mu2e code, but it is required by ROOT persistency
    ExtMonFNALRawHit() : pixelId_(), clock_(), tot_() {}

    const ExtMonFNALPixelId& pixelId() const { return pixelId_; }
    int clock() const { return clock_; }
    unsigned tot() const { return tot_; }

  private:
    ExtMonFNALPixelId pixelId_;
    int clock_;
    unsigned tot_;
  };

  std::ostream& operator<<(std::ostream& os, const ExtMonFNALRawHit& hit);

} // namespace mu2e

#endif /* RecoDataProducts_ExtMonFNALRawHit_hh */
