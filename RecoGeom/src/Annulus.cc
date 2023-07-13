#include "Offline/RecoGeom/inc/Annulus.hh"
#include "Math/VectorUtil.h"
using namespace ROOT::Math::VectorUtil;
namespace mu2e {
  namespace RecoGeom {
    bool Annulus::inBounds(XYZVectorD const& point, double tol) const {
      auto pvec = PerpVector(point - center(),normal());
      return onAnnulus(pvec.R());
    }
  }
}

std::ostream& operator <<(std::ostream& ost, mu2e::RecoGeom::Annulus const& ann) {
  mu2e::RecoGeom::Plane const& plane = static_cast<mu2e::RecoGeom::Plane const&>(ann);
  ost << "Annulus in " << plane << " Inner, Outer radii " << ann.innerRadius() << " , " << ann.outerRadius();
  return ost;
}
