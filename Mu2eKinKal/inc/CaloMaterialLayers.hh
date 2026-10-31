//
// Uniform view of the calorimeter front panel material: the fixed-thickness layers and the cooling pipes, whose
// thickness depends on where the track crosses the pipe plane, and of the crystal case inner support rings. Used both
// to add this material to fits that include a CaloHit and to extrapolate fits without one through the calorimeter.
//
#ifndef Mu2eKinKal_CaloMaterialLayers_hh
#define Mu2eKinKal_CaloMaterialLayers_hh
#include "Offline/KinKalGeom/inc/KinKalGeom.hh"
#include "Offline/KinKalGeom/inc/CaloFrontPanelPipes.hh"
#include "Offline/DataProducts/inc/SurfaceId.hh"
#include "KinKal/Geometry/Plane.hh"
#include "KinKal/Geometry/ParticleTrajectoryIntersect.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <vector>
namespace mu2e {
  struct CaloMaterialLayer {
    std::shared_ptr<KinKal::Plane> surface_;
    SurfaceId sid_;
    std::string const* material_;
    double thickness_; // fixed thickness, for layers
    KKGeom::CaloFrontPanelPipes const* pipes_ = nullptr; // position-dependent thickness, for pipes
    // material depth along the plane normal at a crossing position; 0 means no material is crossed
    double thickness(KinKal::VEC3 const& pos) const { return pipes_ ? pipes_->wallDepth(pos) : thickness_; }
  };

  inline std::vector<CaloMaterialLayer> caloMaterialLayers(KinKalGeom const& kkg) {
    std::vector<CaloMaterialLayer> layers;
    for(auto const& plane : kkg.caloMaterialPlanes())
      layers.push_back(CaloMaterialLayer{plane.surface_,plane.sid_,&plane.material_,plane.thickness_});
    for(auto const& pipes : kkg.caloFrontPanelPipes())
      layers.push_back(CaloMaterialLayer{pipes.plane(),pipes.surfaceId(),&pipes.material(),0.0,&pipes});
    return layers;
  }

  // A crossing of a calorimeter crystal case inner support material, with the path length through that material
  struct CaloSupportXing {
    KKGeom::CaloInnerSupport::Crossing xing_;
    KinKal::Intersection inter_;
    double path_;
    // effective thickness that gives a shell crossing at this intersection the path length through the material
    double thickness() const { return path_*std::fabs(inter_.norm_.Dot(inter_.pdir_)); }
  };

  // Find the support material crossed by a track entering the crystal case. The track must reach the case front moving
  // downstream: if it does so at the crystals it showers there and never crosses the rings; if it does so in the disk
  // hole, it enters the case where it first crosses the inner cylinder. The aluminium and carbon fiber rings span the
  // case length side by side in radius, so the first material is placed where the track enters the case and the second
  // where it then crosses the carbon fiber ring inner radius (outwards into the carbon fiber, or inwards into the
  // aluminium when entering the case front on the carbon fiber), so they are at different times.
  template <class KTRAJ> std::vector<CaloSupportXing> caloSupportXings(KinKal::ParticleTrajectory<KTRAJ> const& ptraj,
      KKGeom::CaloInnerSupport const& support, KinKal::TimeRange const& trange, double tol) {
    using Support = KKGeom::CaloInnerSupport;
    std::vector<CaloSupportXing> xings;
    auto finter = KinKal::intersect(ptraj,*support.frontPlane(),trange,tol);
    if(!finter.good() || finter.pdir_.Z() <= 0.0) return xings;
    double rfront = support.transverseRadius(finter.pos_);
    if(rfront >= support.rcrystal()) return xings;
    bool fromhole = rfront < support.rin();
    auto entry = fromhole ? KinKal::intersect(ptraj,*support.inner(),KinKal::TimeRange(finter.time_,trange.end()),tol)
      : KinKal::intersect(ptraj,*support.front(),trange,tol);
    if(!entry.good()) return xings;
    auto paths = support.pathLengths(entry.pos_,entry.pdir_);
    bool oncf = !fromhole && rfront >= support.rcf();
    auto first = oncf ? Support::carbon : Support::aluminium;
    auto second = oncf ? Support::aluminium : Support::carbon;
    if(paths[first] > 0.0)
      xings.push_back(CaloSupportXing{oncf ? Support::cfFront : fromhole ? Support::alHole : Support::alFront,entry,paths[first]});
    if(paths[second] > 0.0){
      auto cfinter = KinKal::intersect(ptraj,*support.cfInner(),KinKal::TimeRange(entry.time_,trange.end()),tol);
      if(cfinter.good() && cfinter.time_ > entry.time_)
        xings.push_back(CaloSupportXing{oncf ? Support::alRing : Support::cfRing,cfinter,paths[second]});
    }
    return xings;
  }
}
#endif
