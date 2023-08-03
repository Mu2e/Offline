#include "Offline/KinKalGeom/inc/SurfaceMap.hh"
#include "cetlib_except/exception.h"
namespace mu2e {
  SurfaceMap::SurfaceMap() {
    using KinKal::Surface;
    // create the entries for the map; first tracker
    map_[SurfaceId(SurfaceIdEnum::TT_Front)] = std::static_pointer_cast<Surface>(tracker_.frontPtr());
    map_[SurfaceId(SurfaceIdEnum::TT_Mid)] = std::static_pointer_cast<Surface>(tracker_.middlePtr());
    map_[SurfaceId(SurfaceIdEnum::TT_Back)] = std::static_pointer_cast<Surface>(tracker_.backPtr());
    map_[SurfaceId(SurfaceIdEnum::TT_Inner)] = std::static_pointer_cast<Surface>(tracker_.innerPtr());
    map_[SurfaceId(SurfaceIdEnum::TT_Outer)] = std::static_pointer_cast<Surface>(tracker_.outerPtr());
    // DS
    map_[SurfaceId(SurfaceIdEnum::DS_Front)] = std::static_pointer_cast<Surface>(ds_.frontPtr());
    map_[SurfaceId(SurfaceIdEnum::DS_Back)] = std::static_pointer_cast<Surface>(ds_.backPtr());
    map_[SurfaceId(SurfaceIdEnum::DS_Inner)] = std::static_pointer_cast<Surface>(ds_.innerPtr());
    map_[SurfaceId(SurfaceIdEnum::DS_Outer)] = std::static_pointer_cast<Surface>(ds_.outerPtr());
    map_[SurfaceId(SurfaceIdEnum::IPA)] = std::static_pointer_cast<Surface>(ds_.innerProtonAbsorberPtr());
    map_[SurfaceId(SurfaceIdEnum::OPA)] = std::static_pointer_cast<Surface>(ds_.outerProtonAbsorberPtr());
    map_[SurfaceId(SurfaceIdEnum::TSDA)] = std::static_pointer_cast<Surface>(ds_.upstreamAbsorberPtr());
    // target
    map_[SurfaceId(SurfaceIdEnum::ST_Front)] = std::static_pointer_cast<Surface>(st_.frontPtr());
    map_[SurfaceId(SurfaceIdEnum::ST_Back)] = std::static_pointer_cast<Surface>(st_.backPtr());
    map_[SurfaceId(SurfaceIdEnum::ST_Inner)] = std::static_pointer_cast<Surface>(st_.innerPtr());
    map_[SurfaceId(SurfaceIdEnum::ST_Outer)] = std::static_pointer_cast<Surface>(st_.outerPtr());
    for(size_t ifoil=0;ifoil < st_.foils().size();++ifoil){
      map_[SurfaceId(SurfaceIdEnum::ST_Foil,ifoil)] = std::static_pointer_cast<Surface>(st_.foilPtr(ifoil));
    }
    // test CRV; Planes are numbered by their vertical (y) position
    map_[SurfaceId(SurfaceIdEnum::TCRV_Plane,0)] = std::static_pointer_cast<Surface>(tcrv_.t1Ptr());
    map_[SurfaceId(SurfaceIdEnum::TCRV_Plane,1)] = std::static_pointer_cast<Surface>(tcrv_.ex1Ptr());
    map_[SurfaceId(SurfaceIdEnum::TCRV_Plane,2)] = std::static_pointer_cast<Surface>(tcrv_.t2Ptr());

  }
  void SurfaceMap::surfaces(std::vector<SurfaceId> const& ids, std::vector<SurfacePair>& surfs) const {
    surfs.clear();
    surfs.reserve(ids.size());
    for(auto const& id : ids ) {
      auto ifnd = map_.find(id);
      if(ifnd != map_.end()) {
        surfs.push_back(*ifnd);
      } else {
        throw cet::exception("Reco")<<"No Matching Surface"<< std::endl;
      }
    }
  }
}
