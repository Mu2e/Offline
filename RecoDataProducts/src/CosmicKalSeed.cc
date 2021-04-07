#include "RecoDataProducts/inc/CosmicKalSeed.hh"
#include <limits>
namespace mu2e {
  

  HitT0 CosmicKalSeed::t0() const {
    if(segments().size() <= 0) throw cet::exception("RECO")<<"mu2e::CosmicKalSeed: no segments" << std::endl;
// find the segment nearest z=0.  If there's one segment, we are done
    auto iseg = segments().begin();
    if(_segments.size() > 1){
      auto pvel = iseg->state().velocity();
      double vz = pvel.Z(); // to a good approximation, B is along Z
      auto pref = iseg->position3();
      double zmin = pref.Z() + vz*(iseg->tmin()-iseg->tref());
      double zmax = pref.Z() + vz*(iseg->tmax()-iseg->tref());
      if(zmin > 0.0 || zmax < 0.0){
        double mindz = std::min(fabs(zmin), fabs(zmax));
        // find the segment closest to z=0
        auto jseg = ++iseg;
        while(jseg != segments().end()) {
          pvel = jseg->state().velocity();
          vz = pvel.Z();
          pref = jseg->position3();
          zmin = pref.Z() + vz*(jseg->tmin()-jseg->tref());
          zmax = pref.Z() + vz*(jseg->tmax()-jseg->tref());
          double dz = std::min(fabs(zmin),fabs(zmax));
          if(zmin < 0.0 && zmax > 0.0){
            iseg = jseg;
            break;
          } else if(dz < mindz){
            iseg = jseg;
            mindz = dz;
          }
          ++jseg;
        }
      }
    }
    return iseg->t0();
  }
}

