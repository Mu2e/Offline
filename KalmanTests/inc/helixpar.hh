#ifndef helixpar_HH
#define helixpar_HH
#include "Rtypes.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Vector.h"
namespace mu2e
{
  struct helixpar {
    Float_t _d0, _p0, _om, _z0, _td;
    helixpar() : _d0(0.0),_p0(0.0),_om(0.0),_z0(0.0),_td(0.0) {}
    helixpar(const CLHEP::HepVector& pvec) : _d0(pvec[0]),_p0(pvec[1]),_om(pvec[2]),_z0(pvec[3]),_td(pvec[4]) {}
    helixpar(const CLHEP::HepSymMatrix& pcov) : _d0(sqrt(pcov.fast(1,1))),_p0(sqrt(pcov.fast(2,2))),_om(sqrt(pcov.fast(3,3))),
    _z0(sqrt(pcov.fast(4,4))),_td(sqrt(pcov.fast(5,5))) {}
    void reset() { _d0 = _p0 = _om = _z0 = _td = 0.0; }
  };
}
#endif
