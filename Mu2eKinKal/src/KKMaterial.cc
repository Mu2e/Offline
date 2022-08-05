#include "Offline/Mu2eKinKal/inc/KKMaterial.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"

namespace mu2e {
  using MatDBInfo = MatEnv::MatDBInfo;

  KKMaterial::KKMaterial(KKMaterial::Config const& matconfig) :
    filefinder_(matconfig.elements(),matconfig.isotopes(),matconfig.materials()),
    wallmatname_(matconfig.strawWallMaterialName()),
    gasmatname_(matconfig.strawGasMaterialName()),
    wirematname_(matconfig.strawWireMaterialName()),
    eloss_((MatEnv::DetMaterial::energylossmode)matconfig.eloss()),
    matdbinfo_(0) {}

  KKStrawMaterial const& KKMaterial::strawMaterial() const {
    if(matdbinfo_ == 0){
      matdbinfo_ = new MatDBInfo(filefinder_,eloss_);
      Tracker const & tracker = *(GeomHandle<Tracker>());
      auto const& sprop = tracker.strawProperties();
      smat_ = std::make_unique<KKStrawMaterial>(
          sprop,
          matdbinfo_->findDetMaterial(wallmatname_),
          matdbinfo_->findDetMaterial(gasmatname_),
          matdbinfo_->findDetMaterial(wirematname_));
    }
    return *smat_;
  }
}
