#include "Offline/RecoGeom/inc/Ray.hh"
std::ostream& operator <<(std::ostream& ost, mu2e::RecoGeom::Ray const& ray){
  ost << "Ray starting " << ray.start_ << " direction " << ray.dir_;
  return ost;
}
