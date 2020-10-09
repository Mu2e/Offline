#ifndef GeomPrimitives_TubsParams_hh
#define GeomPrimitives_TubsParams_hh

//
// The parameters of a TUBS
//
//
// Original author Rob Kutschke
//

#include <cmath>
#include <ostream>
#include "boost/array.hpp"

#include "CLHEP/Units/PhysicalConstants.h"

namespace mu2e {

  class TubsParams{

  public:

    TubsParams():
      data_()
    {
      data_[0] = 0.;
      data_[1] = 0.;
      data_[2] = 0.;
      data_[3] = 0.;
      data_[4] = 0.;
    }

    TubsParams( double innerRadius,
                double outerRadius,
                double zHalfLength,
                double phi0   = 0.,
                double phiMax = CLHEP::twopi): // spanning angle
      data_()
    {
      data_[0] = innerRadius;
      data_[1] = outerRadius;
      data_[2] = zHalfLength;
      data_[3] = phi0;
      data_[4] = phiMax;
    }

    // Use compiler-generated copy c'tor, copy assignment, and d'tor.

    double innerRadius()  const { return data_[0]; }
    double outerRadius()  const { return data_[1]; }
    double zHalfLength()  const { return data_[2]; }
    double phi0()         const { return data_[3]; }
    double phiMax()       const { return data_[4]; }
    double phiTotal()     const { return data_[4]; }

    double const * data() const { return data_.data(); }

  private:

    boost::array<double,5> data_;

  };

  inline std::ostream& operator<<(std::ostream& ost,
                                  const TubsParams& tp ){
    ost << "("
        << tp.innerRadius() << " "
        << tp.outerRadius() << " "
        << tp.zHalfLength() << " "
        << tp.phi0()        << " "
        << tp.phiMax()      << ")"
        << " )";
    return ost;
  }

}


#endif /* GeomPrimitives_TubsParams_hh */
