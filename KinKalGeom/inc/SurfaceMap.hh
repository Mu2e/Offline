//
//  Factory returning a surface by Id.
//  original author: David Brown (LBNL) 2023
//
#ifndef KinKalGeom_SurfaceMap_hh
#define KinKalGeom_SurfaceMap_hh
#include "Offline/KinKalGeom/inc/Tracker.hh"
#include "Offline/KinKalGeom/inc/StoppingTarget.hh"
#include "Offline/KinKalGeom/inc/DetectorSolenoid.hh"
#include "Offline/KinKalGeom/inc/TestCRV.hh"
#include "Offline/DataProducts/inc/SurfaceId.hh"
#include "KinKal/Geometry/Surface.hh"
#include <memory>
#include <map>
#include <vector>
namespace mu2e {
  class SurfaceMap {
    public:
      using SurfacePtr = std::shared_ptr<KinKal::Surface>;
      using SurfacePair =std::pair<const SurfaceId, SurfacePtr >;
      using SurfacePairCollection = std::vector<SurfacePair>;
      using SurfacePairIter = std::multimap<SurfaceId,SurfacePtr>::const_iterator;
     // default constructor with nominal geometry.  Eventually serve this from GeometryService.  TODO
      SurfaceMap();
      auto const& map() const { return map_; }
      // find a surface by its Id.  Return value is an iterator, which may be null.  Note that if
      // a generic index (-1) is given for surfaces with >1 value, this will return a valid but unspecified matching surface
      auto surface(SurfaceId const& sid) const { return map_.find(sid); }
      // fill a vector of surfaces given a vector of Ids.  This will throw if the id doesn't match
      void surfaces(std::vector<SurfaceId> const& ids, std::vector<SurfacePair>& surfs) const;
      auto const& DS() const {return ds_; }
      auto const& ST() const {return st_; }
      auto const& tracker() const {return tracker_; }
      auto const& TCRV() const {return tcrv_; }
    private:
      // local copy of detector objects; these hold the actual (typed) surface objects
      KinKalGeom::Tracker tracker_;
      KinKalGeom::StoppingTarget st_;
      KinKalGeom::DetectorSolenoid ds_;
      KinKalGeom::TestCRV tcrv_;
      // the actual map
      std::multimap<SurfaceId,SurfacePtr> map_;
  };
}

#endif
