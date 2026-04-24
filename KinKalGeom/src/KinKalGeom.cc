#include "Offline/KinKalGeom/inc/KinKalGeom.hh"
#include "cetlib_except/exception.h"
namespace mu2e {
  KinKalGeom::KinKalGeom() : ProditionsEntity("KinKalGeom") {}

  void KinKalGeom::surfaces(SurfaceId const& id, SurfacePairCollection& surfs) const {
    // find all surfaces that match the input, including index wildcard
    auto frange = map_.equal_range(id);
    for(auto ifnd = frange.first; ifnd != frange.second; ++ifnd){
      surfs.push_back(*ifnd);
    }
  }

  void KinKalGeom::surfaces(SurfaceIdCollection const& ids,SurfacePairCollection& surfs) const {
    surfs.clear();
    surfs.reserve(ids.size());
    for(auto const& id : ids )surfaces(id,surfs);
  }
}
