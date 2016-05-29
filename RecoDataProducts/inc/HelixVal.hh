#ifndef RecoDataProducts_HelixVal_hh
#define RecoDataProducts_HelixVal_hh
//
// helix parameters container
//
// $Id: $
// $Author:  $
// $Date:  $
//
// Original author G. Tassielli
//

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
	double helicity() const { return copysign(1.0,_tanDip*_omega); }

};

} // namespace mu2e

#endif /* RecoDataProducts_HelixVal_hh */
