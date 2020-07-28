#include "TrackerConditions/inc/StrawAlignment.hh"
#include "cetlib_except/exception.h"
#include <cmath>
#include <algorithm>

using namespace std;

namespace mu2e {
  typedef ROOT::Math::XYZVector Vec3; // spatial only vector

  Vec3 StrawAlignment::wireCenter(Straw const& straw, double upos, double hv, Vec3 gravity) const {
  // find the local position
    Vec3 uvw = wireCenterUVW(straw,upos,hv,gravity);
    // translate to global coordintes
    auto midp = straw.getMidPoint();
    auto udir = straw.getDirection(); // U is along the straw
    static CLHEP::Hep3Vector wdir(0.0,0.0,1.0); // W is along Z.  This ignores 2nd order effects due to panel rotation
    auto vdir = wdir.cross(udir);
    midp += udir*uvw.X() + vdir*uvw.Y() + wdir*uvw.Z();
    return  Vec3(midp.x(),midp.y(),midp.z());
  }

  Vec3 StrawAlignment::strawCenter(Straw const& straw, double upos, double hv, Vec3 gravity) const {
  // find the local position
    Vec3 uvw = strawCenterUVW(straw,upos,hv,gravity);
    // translate to global coordintes
    auto midp = straw.getMidPoint();
    auto udir = straw.getDirection(); // U is along the straw
    static CLHEP::Hep3Vector wdir(0.0,0.0,1.0); // W is along Z.  This ignores 2nd order effects due to panel alignment
    auto vdir = wdir.cross(udir);
    midp += udir*uvw.X() + vdir*uvw.Y() + wdir*uvw.Z();
    return  Vec3(midp.x(),midp.y(),midp.z());
  }

  Vec3 StrawAlignment::wireCenterUVW(Straw const& straw, double upos, double hv, Vec3 gravity) const {
    // ignore HV and gravitational sag for now TODO!
    // lookup the wire end data for this straw
    auto iSEA = wireEnds_.find(straw.id());
    if(iSEA == wireEnds_.end())
      throw cet::exception("Reco") << "mu2e::StrawAlignment: no wire end data for straw" << straw.id() << endl;
    return centerUVW(iSEA->second,straw,upos,hv,gravity);
  }

  Vec3 StrawAlignment::strawCenterUVW(Straw const& straw, double upos, double hv, Vec3 gravity) const {
    // ignore HV and gravitational sag for now TODO!
    // lookup the wire end data for this straw
    auto iSEA = strawEnds_.find(straw.id());
    if(iSEA == strawEnds_.end())
      throw cet::exception("Reco") << "mu2e::StrawAlignment: no straw end data for straw" << straw.id() << endl;
    return centerUVW(iSEA->second,straw,upos,hv,gravity);
  }

  Vec3 StrawAlignment::centerUVW(StrawEndAlignment const& sea, Straw const& straw, double upos, double hv, Vec3 gravity) const {
  // U origin is at the straw center.
  // truncate to the length of the physical straw
    upos = std::min(std::max(upos,-straw.halfLength()),straw.halfLength());
    double ulen = 2.0*straw.halfLength();
    double wslope = (sea.dW(StrawEnd(StrawEnd::hv))-sea.dW(StrawEnd(StrawEnd::cal)))/ulen;
    double vslope = (sea.dV(StrawEnd(StrawEnd::hv))-sea.dV(StrawEnd(StrawEnd::cal)))/ulen;
    double wpos = sea.dW(StrawEnd(StrawEnd::cal)) + wslope*(upos+straw.halfLength());
    double vpos = sea.dV(StrawEnd(StrawEnd::cal)) + vslope*(upos+straw.halfLength());
    return Vec3(upos,vpos,wpos);
  }


}


