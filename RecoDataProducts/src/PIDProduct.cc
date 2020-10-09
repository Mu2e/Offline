//
//
//  Original author Vadim Rusu
//

#include "RecoDataProducts/inc/PIDProduct.hh"

using namespace std;

namespace mu2e {

  // Constructors
  PIDProduct::PIDProduct() {}


  PIDProduct::PIDProduct(const PIDProduct & p) {

    _residualsSlope = p._residualsSlope;
    _residualsSlopeError = p._residualsSlopeError;
    _trkid = p._trkid;
    _logeprob = p._logeprob;
    _logmprob = p._logmprob;

  }

  // operator overloading
  PIDProduct & PIDProduct::operator= (const PIDProduct & p) {
    _residualsSlope = p._residualsSlope;
    _residualsSlopeError = p._residualsSlopeError;
    _logeprob = p._logeprob;
    _logmprob = p._logmprob;
    _trkid = p._trkid;
    return (*this);
  }


  void PIDProduct::clear () {
    _residualsSlope = -999.;
    _residualsSlopeError = -999.;
    _trkid = -999.;
    _logeprob = -999.;
    _logmprob = -999.;

  }




} // end namespace mu2e


