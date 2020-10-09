#ifndef TrackerConditions_StrawEndAlignment_hh
#define TrackerConditions_StrawEndAlignment_hh
//
// Class to specify the alignment of the ends of an individual straw or wire WRT the nominal geometry.
// This is currently just a displacement: angles at the end might eventually be needed, at least for straws.
// End alignment is necessary but insufficient to describe the straw or wire as an object in space, as
// gravitational and electrostatic forces act on the straw and wire and alter its position along the length.  The proditions
// object StrawAlignment takes this data plus tension and gravity and electrostatics, and computes
// the full shape of the straw or wire in space
//
// Note this object uses the straw-local coordinate system UVW, where U points along the straw (from Cal to HV),
// W points along the nominal Mu2e Z direction, and V is perpendicular to those
// 
// Original author David Brown (7/2020)
//
// Mu2e includes
#include "DataProducts/inc/StrawId.hh"
#include "DataProducts/inc/StrawEnd.hh"
#include "Math/SMatrix.h"
namespace mu2e {
  struct StrawEndAlignment {  
  typedef ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepSym<double,2> > CMAT;
    StrawId id_; // defines which straw.
    float dV_[2]; // displacement in V direction at each end
    float dW_[2]; // displacement in W direction at each end
    CMAT covar_; // covariance matrix of measurements of displacements
    // interface for alignment
    float dV(StrawEnd end) const { return dV_[end.end()]; }
    float dW(StrawEnd end) const { return dW_[end.end()]; }
  };
}
#endif
