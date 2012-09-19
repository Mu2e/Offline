// $Id: ExtMonFNALTrkParam.hh,v 1.1 2012/09/19 03:54:19 gandr Exp $
// $Author: gandr $
// $Date: 2012/09/19 03:54:19 $
//
// Original author Andrei Gaponenko
//

#ifndef RecoDataProducts_ExtMonFNALTrkParam_hh
#define RecoDataProducts_ExtMonFNALTrkParam_hh

#include <ostream>

//#define MATRIX_BOUND_CHECK /* for CLHEP matrix classes */
#include "CLHEP/Matrix/Vector.h"
//#undef MATRIX_BOUND_CHECK

namespace mu2e {

  //================================================================
  class ExtMonFNALTrkParam {
  public:
    enum ParIndex {POSX, SLOPEX, POSY, SLOPEY, RINV, NPARS};

    // The Z axis has a kink at the detector reference point z==0.
    // It follows the axis of the upstream stack for z>0,
    // and of the downstream stack for z<0.

    // A track is a straight line in each of the stacks.
    // The section in the magnet is a piece of helix.
    // RINV == 1/Rtrack in projection on the YZ plane.
    // We define a track by its straight line parameters at a given
    // z=z0 and 1/R.

    double z0() const { return z0_; }

    // Parameters of the straight line equation
    //   x = posx + slopex*(z-z0).
    double posx() const { return pars_[POSX]; }
    double slopex() const { return pars_[SLOPEX]; }

    // Parameters in the YZ view
    //   y = posy + slopey*(z-z0).
    double posy() const { return pars_[POSY]; }
    double slopey() const { return pars_[SLOPEY]; }

    // 1/R
    double rinv() const { return pars_[RINV]; }

    const CLHEP::HepVector& pars() const { return pars_; }

    ExtMonFNALTrkParam() : pars_(NPARS, 0.), z0_() {}
    ExtMonFNALTrkParam(double z0, const CLHEP::HepVector& pars)
      : pars_(pars), z0_(z0)
    {}

    void setz0(double z) { z0_ = z; }
    void setposx(double x) { pars_[POSX] = x; }
    void setslopex(double sx) { pars_[SLOPEX] = sx; }
    void setposy(double y) { pars_[POSY] = y; }
    void setslopey(double sy) { pars_[SLOPEY] = sy; }
    void setrinv(double rinv) { pars_[RINV] = rinv; }

  private:
    CLHEP::HepVector pars_;
    double z0_;
  };

  std::ostream& operator<<(std::ostream& os, const ExtMonFNALTrkParam& c);

} // namespace mu2e

#endif /* RecoDataProducts_ExtMonFNALTrkParam_hh */
