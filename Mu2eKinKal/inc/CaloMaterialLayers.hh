//
// Uniform view of the calorimeter front panel material: the fixed-thickness layers and the cooling pipes, whose
// thickness depends on where the track crosses the pipe plane. Used both to add the front panel material to fits
// that include a CaloHit and to extrapolate fits without one through the front panels.
//
#ifndef Mu2eKinKal_CaloMaterialLayers_hh
#define Mu2eKinKal_CaloMaterialLayers_hh
#include "Offline/KinKalGeom/inc/KinKalGeom.hh"
#include "Offline/KinKalGeom/inc/CaloFrontPanelPipes.hh"
#include "Offline/DataProducts/inc/SurfaceId.hh"
#include "KinKal/Geometry/Plane.hh"
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
}
#endif
