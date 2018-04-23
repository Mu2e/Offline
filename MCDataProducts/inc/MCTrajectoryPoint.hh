#ifndef MCDataProducts_MCTrajectoryPoint_hh
#define MCDataProducts_MCTrajectoryPoint_hh

#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class MCTrajectoryPoint {
  private:
    float x_;
    float y_;
    float z_;
    float t_;
    float kineticEnergy_;

  public:
    // default ctr required for ROOT I/O
    MCTrajectoryPoint(): x_(), y_(), z_(), t_(), kineticEnergy_() {}

    MCTrajectoryPoint(const CLHEP::Hep3Vector& pos, float t, float kE)
      : x_(pos.x()), y_(pos.y()), z_(pos.z()), t_(t), kineticEnergy_(kE)
    {}

    CLHEP::Hep3Vector pos() const { return CLHEP::Hep3Vector(x_, y_, z_); }

    float x() const { return x_; }
    float y() const { return y_; }
    float z() const { return z_; }
    float t() const { return t_; }
    float kineticEnergy() const { return kineticEnergy_; }
  };
}

#endif /*MCDataProducts_MCTrajectoryPoint_hh*/
