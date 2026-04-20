//
//  Top-level object for referencing geometry (surfaces) used in KinKal, either
//  hierarchically or by Surface Id
//  original author: David Brown (LBNL) 2023
//
#ifndef KinKalGeom_KinKalGeom_hh
#define KinKalGeom_KinKalGeom_hh
#include "Offline/KinKalGeom/inc/Tracker.hh"
#include "Offline/KinKalGeom/inc/StoppingTarget.hh"
#include "Offline/KinKalGeom/inc/DetectorSolenoid.hh"
#include "Offline/KinKalGeom/inc/CRV.hh"
#include "Offline/KinKalGeom/inc/TestCRV.hh"
#include "Offline/DataProducts/inc/SurfaceId.hh"
#include "KinKal/Geometry/Surface.hh"
#include <memory>
#include <map>
#include <vector>
namespace mu2e {
  class KinKalGeom {
    public:
      using SurfacePtr = std::shared_ptr<KinKal::Surface>;
      using SurfacePair =std::pair<const SurfaceId, SurfacePtr >;
      using SurfacePairCollection = std::vector<SurfacePair>;
      using SurfacePairIter = std::multimap<SurfaceId,SurfacePtr>::const_iterator;
      using KKGMap = std::multimap<SurfaceId,SurfacePtr>;
      // default constructor, now using GeometryService. This
      KinKalGeom(){}
      auto const& map() const { check_init(); return map_; }
      // find a surface by its Id.  Return value is an iterator, which may be null.  Note that if
      // a generic index (-1) is given for surfaces with >1 value, this will return a valid but unspecified matching surface
      auto surface(SurfaceId const& sid) const { return map_.find(sid); }
      // fill a vector of surfaces given a vector of Ids.  This will throw if an id doesn't match
      void surfaces(std::vector<SurfaceId> const& ids, std::vector<SurfacePair>& surfs) const;
      // hierarchical accessors
      auto const& DS() const {return ds_; }
      auto const& ST() const {return st_; }
      auto const& tracker() const {return tracker_; }
      auto const& CRV() const {return crv_; }
      auto const& TCRV() const {return tcrv_; }
    private:
      mutable bool initialized_ = false; // defer construction to allow services to be established
      void check_init() const;
      void initialize();
      // local copy of detector objects; these hold the actual (typed) surface objects
      KKGeom::Tracker tracker_;
      KKGeom::StoppingTarget st_;
      KKGeom::DetectorSolenoid ds_;
      KKGeom::CRV crv_;
      KKGeom::TestCRV tcrv_;
      // the map used to find surfaces by Id
      KKGMap map_;
  };
}

#endif
