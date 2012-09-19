#ifndef RecoDataProducts_ExtMonFNALRawHit_hh
#define RecoDataProducts_ExtMonFNALRawHit_hh
//
// Unpacked data from the ExtMonFNAL pixel detector.
//
// $Id: ExtMonFNALRawHit.hh,v 1.3 2012/09/19 04:59:00 gandr Exp $
// $Author: gandr $
// $Date: 2012/09/19 04:59:00 $
//
// Original author Andrei Gaponenko
//

#include <ostream>
#include <cassert>

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

#ifdef ENABLE_MU2E_GENREFLEX_HACKS
  // FIXME: Hits can't be compared, but the combination of genreflex and
  // current art::PtrVector insists that on defining this operator.
  bool operator<(const ExtMonFNALRawHit&, const ExtMonFNALRawHit&) {
    assert(false);
    return false;
  }
#endif/*ENABLE_MU2E_GENREFLEX_HACKS*/

} // namespace mu2e

#endif /* RecoDataProducts_ExtMonFNALRawHit_hh */
