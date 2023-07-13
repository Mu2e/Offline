#include "Offline/RecoGeom/inc/Disk.hh"

std::ostream& operator <<(std::ostream& ost, mu2e::RecoGeom::Disk const& disk) {
  mu2e::RecoGeom::Plane const& plane = static_cast<mu2e::RecoGeom::Plane const&>(disk);
  ost << "Disk in " << plane << " radius " << disk.outerRadius();
  return ost;
}
