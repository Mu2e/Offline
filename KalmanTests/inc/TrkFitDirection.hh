// class to define mu2e track fit direction convention.  In a Kalman fit, this is the direction the partlcle
// physically travels WRT time, and the pitch sign.  Together with the BField and the particle charge,
// this also defines the angular velocity sign.
//
// $Id: TrkFitDirection.hh,v 1.1 2012/07/25 20:56:57 brownd Exp $
// $Author: brownd $ 
// $Date: 2012/07/25 20:56:57 $
#ifndef TrkFitDirection_HH
#define TrkFitDirection_HH
#include <string>

namespace mu2e 
{
  class TrkFitDirection {
    public:
// define the fit direction as downstream (towards positive Z) or upstream (towards negative Z).
      enum FitDirection {downstream=0,upstream};
      TrkFitDirection(FitDirection fdir=downstream);
      FitDirection fitDirection() const { return _fdir; }
      std::string const& name() const;
    private:
      FitDirection _fdir;
  };
}
#endif
