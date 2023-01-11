#include "Offline/TrackerConditions/inc/DriftInfo.hh"
namespace mu2e {
  double DriftInfo::maxdvar_(2.4*2.4/3.0); // effective maximum drift distance variance: should come from StrawResponse TODO

  double DriftInfo::nullDistanceVar() const {
    double rdrift = std::max(0.0,driftDistance_); // don't count negative drift distances
    return rdrift*rdrift + nullDistanceError_*nullDistanceError_;
  }
}
