#ifndef RecoDataProducts_HelixVal_hh
#define RecoDataProducts_HelixVal_hh
//
// BaBar definition helix parameters.  This is a mixed geometric/kinematic helix as the
// signs include time propagation information
//

// Mu2e
#include "DataProducts/inc/Helicity.hh"
// BTrk
#include "BTrk/TrkBase/HelixParams.hh"
// CLHEP
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"
// C includes 
#include <Rtypes.h>
#include <math.h>

namespace mu2e {

  struct HelixVal {
    HelixVal() { for(size_t ipar=0;ipar<5;++ipar)_pars[ipar] = 0;}
    HelixVal(CLHEP::HepVector const& pvec);
    HelixVal& operator = (CLHEP::HepVector const& pvec);

    Float_t d0() const { return _pars[HelixParams::d0Index]; }
    Float_t phi0() const { return _pars[HelixParams::phi0Index]; }
    Float_t omega() const { return _pars[HelixParams::omegaIndex]; }
    Float_t z0() const { return _pars[HelixParams::z0Index]; }
    Float_t tanDip() const { return _pars[HelixParams::tanDipIndex]; }

    Float_t& d0() { return _pars[HelixParams::d0Index]; }
    Float_t& phi0() { return _pars[HelixParams::phi0Index]; }
    Float_t& omega() { return _pars[HelixParams::omegaIndex]; }
    Float_t& z0() { return _pars[HelixParams::z0Index]; }
    Float_t& tanDip() { return _pars[HelixParams::tanDipIndex]; }
// convert to CLHEP.
    void hepVector(CLHEP::HepVector& pvec) const;
// this should really be a HepVector, or whatever algebraic vector class we use, FIXME!
    Float_t _pars[5];
    // helicity is given by the product of the signs of tandip (axial motion) and omega (angular momentum)
    Helicity helicity() const { return Helicity(tanDip()*omega()); }
  };

  struct HelixCov {
// lower-diagonal array.  Unfortunately genreflex can't handle std::Array, so I must
// hard-code the limit
    Float_t _cov[15];
    HelixCov();
    HelixCov(CLHEP::HepSymMatrix const& pcov);
    HelixCov& operator = (CLHEP::HepSymMatrix const& pcov);
    void symMatrix(CLHEP::HepSymMatrix& pcov) const;
  };

} // namespace mu2e

#endif /* RecoDataProducts_HelixVal_hh */
