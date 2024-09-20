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
    stmatname_(matconfig.STMaterialName()),
    eloss_((DetMaterial::energylossmode)matconfig.eloss()),
    matdbinfo_(std::make_unique<MatDBInfo>(filefinder_,eloss_)) {
      // std::cout << "Setting KinKal scattering Dahl-Lynch fraction to " << matconfig.dahlLynchScatteringFraction() << std::endl;
      // MatEnv forces const-cast to change scattering parameter
      const_cast<DetMaterial*>(matdbinfo_->findDetMaterial(wallmatname_))->setScatterFraction(matconfig.dahlLynchScatteringFraction());
      const_cast<DetMaterial*>(matdbinfo_->findDetMaterial(gasmatname_))->setScatterFraction(matconfig.dahlLynchScatteringFraction());
      const_cast<DetMaterial*>(matdbinfo_->findDetMaterial(wirematname_))->setScatterFraction(matconfig.dahlLynchScatteringFraction());
      const_cast<DetMaterial*>(matdbinfo_->findDetMaterial(ipamatname_))->setScatterFraction(matconfig.dahlLynchScatteringFraction());
      const_cast<DetMaterial*>(matdbinfo_->findDetMaterial(stmatname_))->setScatterFraction(matconfig.dahlLynchScatteringFraction());
    }

  KKStrawMaterial const& KKMaterial::strawMaterial() const {
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
