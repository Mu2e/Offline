#ifndef RecoDataProducts_HelixVal_hh
#define RecoDataProducts_HelixVal_hh
//
// BaBar definition helix parameters.  This is a mixed geometric/kinematic helix as the
// signs include time propagation information
//

// Mu2e
#include "DataProducts/inc/Helicity.hh"
// C includes 
#include <math.h>

namespace mu2e {

struct HelixVal {
        HelixVal() :
                _d0    (0),
                _phi0  (0),
                _omega (0),
                _z0    (0),
                _tanDip(0){}

        double                _d0;
        double                _phi0;
        double                _omega;
        double                _z0;
        double                _tanDip;

// helicity is given by the product of the signs of tandip (axial motion) and omega (angular momentum)
	Helicity helicity() const { return Helicity(_tanDip*_omega); }

};

} // namespace mu2e

#endif /* RecoDataProducts_HelixVal_hh */
