//
// BaBar convention helix parameters.  This struct is for us in TTrees, the actual helix fit parameters are
// defined by BaBar classes
// Dave Brown (LBNL)
//
#ifndef helixpar_HH
#define helixpar_HH
#include "Rtypes.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Vector.h"
#include "RecoDataProducts/inc/HelixVal.hh"
namespace mu2e
{
  struct helixpar {
    Float_t _d0; // transverse distance of closest approach to the z axis, signed by the track angular momentum
    Float_t _p0; // aximuth of the momentum vector at the point of closest approach to the z axis
    Float_t _om; // curvature (inverse radius), signed by the charge.  This is position independent
    Float_t _z0; // z position of the track at the transverse point of closest approach
    Float_t _td; // tangent of the dip angle (= pi/2 - theta).  This is independent of position
    helixpar() : _d0(0.0),_p0(0.0),_om(0.0),_z0(0.0),_td(0.0) {}
    helixpar(const CLHEP::HepVector& pvec) : _d0(pvec[0]),_p0(pvec[1]),_om(pvec[2]),_z0(pvec[3]),_td(pvec[4]) {}
    helixpar(const CLHEP::HepSymMatrix& pcov) : _d0(sqrt(pcov.fast(1,1))),_p0(sqrt(pcov.fast(2,2))),_om(sqrt(pcov.fast(3,3))),
    _z0(sqrt(pcov.fast(4,4))),_td(sqrt(pcov.fast(5,5))) {}
    helixpar(HelixVal const& hval ) : _d0(hval.d0()), _p0(hval.phi0()), _om(hval.omega()), _z0(hval.z0()), _td(hval.tanDip()) {}
    void reset() { _d0 = _p0 = _om = _z0 = _td = 0.0; }
    static std::string leafnames() { static std::string leaves; leaves = std::string("d0/F:p0/F:om/F:z0/F:td/F"); return leaves; }
  };
}
#endif
