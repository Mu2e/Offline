//
// BaBar definition helix parameters.  This is a mixed geometric/kinematic helix as the
// signs include time propagation information
//
#include "RecoDataProducts/inc/HelixVal.hh"
using CLHEP::HepVector;
using CLHEP::HepSymMatrix;

namespace mu2e {
  HelixVal::HelixVal(HepVector const& pvec) {
    for(int irow=0; irow < HelixParams::NHLXPRM; ++irow)
      _pars[irow] = pvec[irow];
  }

  HelixVal& HelixVal::operator = (HepVector const& pvec) {
    for(int irow=0; irow < HelixParams::NHLXPRM; ++irow)
      _pars[irow] = pvec[irow];
    return *this;
  }

  void HelixVal::hepVector(HepVector& pvec) const {
    pvec = HepVector(HelixParams::NHLXPRM,0);
    for(int irow=0; irow < HelixParams::NHLXPRM; ++irow)
      pvec[irow] = _pars[irow];
  }

  HelixCov::HelixCov()  {
    for(size_t ipar=0;ipar<15;++ipar)
      _cov[ipar]=0;
  }

  HelixCov::HelixCov(HepSymMatrix const& pcov) {
    size_t index(0);
    for(int irow=0; irow < HelixParams::NHLXPRM; ++irow){
      for(int icol=0; icol <= irow; ++icol){
	_cov[index] = pcov(irow+1,icol+1); // CLHEP convention!
	++index;
      }
    }
  }

  HelixCov& HelixCov::operator = (HepSymMatrix const& pcov) {
    size_t index(0);
    for(int irow=0; irow < HelixParams::NHLXPRM; ++irow){
      for(int icol=0; icol <= irow; ++icol){
	_cov[index] = pcov(irow+1,icol+1); // CLHEP convention!
	++index;
      }
    }
    return *this;
  }

  void HelixCov::symMatrix(HepSymMatrix& pcov) const {
  // dimension the matrix appropriately
    if(pcov.num_row() != 5)pcov = HepSymMatrix(5,0);
    size_t index(0);
    for(int irow=0; irow < HelixParams::NHLXPRM; ++irow){
      for(int icol=0; icol <= irow; ++icol){
	pcov(irow+1,icol+1) =  _cov[index]; // CLHEP convention!
	++index;
      }
    }
  }
  
  // helix geometry functions
  float HelixVal::phi(float fltlen) const {
    return phi0() + omega()*fltlen*cosDip();
  }

  void HelixVal::direction(float fltlen, XYZVec& dir) const {
    float phival = phi(fltlen);
    float cd = cosDip();
    float sd = sinDip();
    dir = XYZVec(cd*cos(phival),cd*sin(phival),sd);
  }

  void HelixVal::position(float fltlen, XYZVec& pos) const {
    float phival = phi(fltlen);
    float invomega = 1.0/omega();
    float rval = invomega + d0(); 
    pos = XYZVec(invomega*sin(phival) - rval*sin(phi0()),
	-invomega*cos(phival) + rval*cos(phi0()),
	z0() + fltlen*sinDip());
  }

  void HelixVal::position(const XYZVec& pos, float& fltlen) const {
    float z = pos.z();
    fltlen = (z - z0()) / sinDip();
  }

}
