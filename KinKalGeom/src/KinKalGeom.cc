#include "Offline/KinKalGeom/inc/KinKalGeom.hh"
#include "cetlib_except/exception.h"
namespace mu2e {

  void KinKalGeom::check_init() const {
    if(!initialized_){
      // cast off const
      (const_cast<KinKalGeom*>(this))->initialize();
      initialized_ = true;
    }
  }

  void KinKalGeom::initialize() {
    using KinKal::Surface;
    // create the entries for the map; first tracker
    tracker_.addSurfaces(map_);
    // DS
    map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::DS_Front),std::static_pointer_cast<Surface>(ds_.frontPtr())));
    map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::DS_Back),std::static_pointer_cast<Surface>(ds_.backPtr())));
    map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::DS_Inner),std::static_pointer_cast<Surface>(ds_.innerPtr())));
    map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::DS_Outer),std::static_pointer_cast<Surface>(ds_.outerPtr())));
    map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::IPA),std::static_pointer_cast<Surface>(ds_.innerProtonAbsorberPtr())));
    map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::IPA_Front),std::static_pointer_cast<Surface>(ds_.innerProtonAbsorberFrontPtr())));
    map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::IPA_Back),std::static_pointer_cast<Surface>(ds_.innerProtonAbsorberBackPtr())));
    map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::OPA),std::static_pointer_cast<Surface>(ds_.outerProtonAbsorberPtr())));
    map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::TSDA),std::static_pointer_cast<Surface>(ds_.upstreamAbsorberPtr())));
    // target
    map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::ST_Front),std::static_pointer_cast<Surface>(st_.frontPtr())));
    map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::ST_Back),std::static_pointer_cast<Surface>(st_.backPtr())));
    map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::ST_Inner),std::static_pointer_cast<Surface>(st_.innerPtr())));
    map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::ST_Outer),std::static_pointer_cast<Surface>(st_.outerPtr())));
    for(size_t ifoil=0;ifoil < st_.foils().size();++ifoil){
      map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::ST_Foils,ifoil),std::static_pointer_cast<Surface>(st_.foilPtr(ifoil))));
    }
    // test CRV; Planes are numbered by their vertical (y) position
    map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::TCRV,0),std::static_pointer_cast<Surface>(tcrv_.t1Ptr())));
    map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::TCRV,1),std::static_pointer_cast<Surface>(tcrv_.ex1Ptr())));
    map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::TCRV,2),std::static_pointer_cast<Surface>(tcrv_.t2Ptr())));
  }
  void KinKalGeom::surfaces(SurfaceIdCollection const& ids,SurfacePairCollection& surfs) const {
    surfs.clear();
    surfs.reserve(ids.size());
    for(auto const& id : ids ) {
      // find all surfaces whose keys match
      auto frange = map_.equal_range(id);
      for(auto ifnd = frange.first; ifnd != frange.second; ++ifnd){
        surfs.push_back(*ifnd);
      }
    }
  }
}
