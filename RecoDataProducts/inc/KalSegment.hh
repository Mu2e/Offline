//
//  Class representing a segment of the Kalman filter track fit.  The segment contains part of the fit
//  result, accurate and unchanging over the specified range of the particle trajectory.
//  Segments can be strung together to represent the full Kalman fit, or they can be sampled
//  to provide valid fit results in subset of positions.
//  Original Author: Dave Brown (LBNL) 31 Aug. 2016
//
#ifndef RecoDataProducts_KalSegment_HH
#define RecoDataProducts_KalSegment_HH
#include <Rtypes.h>
#include "RecoDataProducts/inc/HelixVal.hh"
#include "DataProducts/inc/XYZVec.hh"

namespace mu2e {
  struct KalSegment {
    KalSegment() : _fmin(0.0), _fmax(-1.0), _dflt(0.0), _mom(-1.0), _momerr(-1.0) {}
    Float_t fmin() const { return _fmin; } // local helix flight limits
    Float_t fmax() const { return _fmax; }
    HelixVal const& helix() const { return _helix; }
    HelixCov const& covar() const { return _hcov; }
    Float_t mom() const { return _mom; }
    Float_t momerr() const { return _momerr; }
    void mom(float fltlen, XYZVec& momvec) const { helix().direction(fltlen,momvec); momvec *= mom(); } // momentum as a function of local flightlength
    // translate back and forth to global and local flightlengths
    float localFlt(float globalflt) const { return globalflt + _dflt; }
    float globalFlt(float localflt) const { return localflt - _dflt; }
    Float_t _fmin, _fmax; // local Flight lengths along the helix for which this segment is valid.
    Float_t _dflt; // difference between local and global flight length
    HelixVal _helix; // helix valid for this segment
    HelixCov _hcov; // covariance matrix for this helix
    Float_t _mom; // scalar momentum of this helix segment
    Float_t _momerr; // error on the scalar momentum
  };
}
#endif

