#include "Offline/Mu2eKinKal/inc/KKMaterial.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "KinKal/MatEnv/DetMaterial.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"

namespace mu2e {
  using MatDBInfo = MatEnv::MatDBInfo;
  using MatEnv::DetMaterial;

  KKMaterial::KKMaterial(KKMaterial::Config const& matconfig) :
    filefinder_(matconfig.elements(),matconfig.isotopes(),matconfig.materials()),
    wallmatname_(matconfig.strawWallMaterialName()),
    gasmatname_(matconfig.strawGasMaterialName()),
    wirematname_(matconfig.strawWireMaterialName()),
    ipamatname_(matconfig.IPAMaterialName()),
    stmatname_(matconfig.STMaterialName()) {
      MatEnv::DetMaterialConfig dmconf;
      dmconf.elossmode_ = (DetMaterial::energylossmode)matconfig.elossMode();
      dmconf.scatterfrac_solid_ = matconfig.solidScatter();
      dmconf.scatterfrac_gas_ = matconfig.gasScatter();
      dmconf.ebrehmsfrac_ = matconfig.eBrehms();
      matdbinfo_ = std::make_unique<MatDBInfo>(filefinder_,dmconf);
    }

  KKStrawMaterial const& KKMaterial::strawMaterial() const {
    // deferred construction as this object depends on the tracker, which is created at beginJob
    if(smat_ == nullptr){
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
