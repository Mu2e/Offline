// class to define mu2e track fit direction convention.  In a Kalman fit, this is the direction the partlcle
// physically travels WRT time, and the pitch sign.  Together with the BField and the particle charge,
// this also defines the angular velocity sign.
//
#ifndef TrkFitDirection_HH
#define TrkFitDirection_HH
#include <string>
#include <cmath>

namespace mu2e
{
  class TrkFitDirection {
    public:
// define the fit direction as downstream (towards positive Z) or upstream (towards negative Z).
      enum FitDirection {downstream=0,upstream,unknown};
      TrkFitDirection(FitDirection fdir=downstream);
      //accessors
      FitDirection fitDirection() const { return _fdir; }
      // return the SIGN of the z component of velocity (magnitude is not returned)
      double dzdt() const { return _fdir == downstream ? 1.0 : _fdir == upstream ? -1.0 : 0.0; }
      std::string const& name() const;
      bool operator == ( TrkFitDirection const& other) const { return _fdir == other._fdir; }
      bool operator != ( TrkFitDirection const& other) const { return _fdir != other._fdir; }
    private:
      FitDirection _fdir;
  };
}
#endif
