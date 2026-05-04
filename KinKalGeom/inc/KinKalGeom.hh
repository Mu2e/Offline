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
#include "Offline/Mu2eInterfaces/inc/Detector.hh"
#include "Offline/Mu2eInterfaces/inc/ProditionsEntity.hh"
#include <memory>
#include <map>
#include <vector>
namespace mu2e {
  class KinKalGeom : public Detector, public ProditionsEntity {
    public:
      using SurfacePtr = std::shared_ptr<KinKal::Surface>;
      using SurfacePair =std::pair<const SurfaceId, SurfacePtr >;
      using SurfacePairCollection = std::vector<SurfacePair>;
      using SurfacePairIter = std::multimap<SurfaceId,SurfacePtr>::const_iterator;
      using KKGMap = std::multimap<SurfaceId,SurfacePtr>;
      // default constructor, now using GeometryService
      KinKalGeom();
      // accessor to the raw map
      auto const& map() const { return map_; }
      // find all surfaces that match an Id. Return vector can have >1 entry if the wildcard index (-1) is provided
      void surfaces(SurfaceId const& sid, std::vector<SurfacePair>& surfs) const;
      // find all surfaces that match a vector of Ids
      void surfaces(std::vector<SurfaceId> const& ids, std::vector<SurfacePair>& surfs) const;
      // hierarchical accessors
      auto const& DS() const {return ds_; }
      auto const& ST() const {return st_; }
      auto const& tracker() const {return tracker_; }
//      auto const& CRV() const {return crv_; }
      auto const& TCRV() const {return tcrv_; }
    private:
      // local copy of detector objects; these hold the actual (typed) surface objects
      std::unique_ptr<KKGeom::Tracker> tracker_;
      std::unique_ptr<KKGeom::DetectorSolenoid> ds_;
      std::unique_ptr<KKGeom::StoppingTarget> st_;
      //KKGeom::CRV crv_;
      std::unique_ptr<KKGeom::TestCRV> tcrv_;
      // the map used to find surfaces by Id
      KKGMap map_;
      // allow GeometryService to access internals, to avoid link loop
      friend class KinKalGeomMaker;
  };
}

#endif
