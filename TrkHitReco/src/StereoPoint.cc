#include "Offline/TrkHitReco/inc/StereoPoint.hh"

std::ostream& operator << (std::ostream& ost, mu2e::StereoPoint const& pt) {
  ost << "StereoPoint z " << pt.z() << " " << pt.point();
  return ost;
}
