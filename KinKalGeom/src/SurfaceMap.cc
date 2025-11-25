#include "Offline/KinKalGeom/inc/SurfaceMap.hh"
#include "cetlib_except/exception.h"
namespace mu2e {
  SurfaceMap::SurfaceMap() {
    using KinKal::Surface;
    // create the entries for the map; first tracker
    map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::TT_Front),std::static_pointer_cast<Surface>(tracker_.frontPtr())));
    map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::TT_Mid),std::static_pointer_cast<Surface>(tracker_.middlePtr())));
    map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::TT_Back),std::static_pointer_cast<Surface>(tracker_.backPtr())));
    map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::TT_Inner),std::static_pointer_cast<Surface>(tracker_.innerPtr())));
    map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::TT_Outer),std::static_pointer_cast<Surface>(tracker_.outerPtr())));
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

    //calo FIXME - need to remove the frontPtr
 map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::EMC_Disk_0_SurfIn),std::static_pointer_cast<Surface>(calo_.EMC_Disk_0_SurfInPtr()))); map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::EMC_Disk_0_SurfOut),std::static_pointer_cast<Surface>(calo_.EMC_Disk_0_SurfOutPtr()))); map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::EMC_Disk_1_SurfIn),std::static_pointer_cast<Surface>(calo_.EMC_Disk_1_SurfInPtr()))); map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::EMC_Disk_1_SurfOut),std::static_pointer_cast<Surface>(calo_.EMC_Disk_1_SurfOutPtr()))); map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::EMC_Disk_0_EdgeIn),std::static_pointer_cast<Surface>(calo_.EMC_Disk_0_EdgeInPtr()))); map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::EMC_Disk_0_EdgeOut),std::static_pointer_cast<Surface>(calo_.EMC_Disk_0_EdgeOutPtr()))); map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::EMC_Disk_1_EdgeIn),std::static_pointer_cast<Surface>(calo_.EMC_Disk_1_EdgeInPtr()))); map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::EMC_Disk_1_EdgeOut),std::static_pointer_cast<Surface>(calo_.EMC_Disk_1_EdgeOutPtr()))); map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::EMC_0_FrontIn),std::static_pointer_cast<Surface>(calo_.EMC_0_FrontInPtr()))); map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::EMC_0_FrontOut),std::static_pointer_cast<Surface>(calo_.EMC_0_FrontOutPtr()))); map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::EMC_1_FrontIn),std::static_pointer_cast<Surface>(calo_.EMC_1_FrontInPtr()))); map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::EMC_1_FrontOut),std::static_pointer_cast<Surface>(calo_.EMC_1_FrontOutPtr()))); map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::EMC_2_FrontIn),std::static_pointer_cast<Surface>(calo_.EMC_2_FrontInPtr()))); map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::EMC_2_FrontOut),std::static_pointer_cast<Surface>(calo_.EMC_2_FrontOutPtr()))); map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::EMC_3_FrontIn),std::static_pointer_cast<Surface>(calo_.EMC_3_FrontInPtr()))); map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::EMC_3_FrontOut),std::static_pointer_cast<Surface>(calo_.EMC_3_FrontOutPtr())));

  }
  void SurfaceMap::surfaces(SurfaceIdCollection const& ids,SurfacePairCollection& surfs) const {
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
