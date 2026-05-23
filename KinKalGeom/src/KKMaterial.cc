#include "Offline/KinKalGeom/inc/KKMaterial.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "KinKal/MatEnv/DetMaterial.hh"

namespace mu2e {
  using MatDBInfo = MatEnv::MatDBInfo;
  using MatEnv::DetMaterial;

  KKMaterial::KKMaterial(KKMaterial::Config const& matconfig,Tracker const& tracker) :
    elementsBaseName_(matconfig.elements()),
    isotopesBaseName_(matconfig.isotopes()),
    materialsBaseName_(matconfig.materials()),
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
      matdbinfo_ = std::make_unique<MatDBInfo>(*this,dmconf);
      smat_ = std::make_unique<KKStrawMaterial>(
          tracker.strawProperties(),
          matdbinfo_->findDetMaterial(wallmatname_),
          matdbinfo_->findDetMaterial(gasmatname_),
          matdbinfo_->findDetMaterial(wirematname_));
    }

// implement file finder policy using mu2e file policy
  std::string KKMaterial::findFile( std::string const& basename ) const {
    ConfigFileLookupPolicy policy;
    return policy( basename ); }
  std::string KKMaterial::matElmDictionaryFileName() const { return findFile(elementsBaseName_); }
  std::string KKMaterial::matIsoDictionaryFileName() const { return findFile(isotopesBaseName_); }
  std::string KKMaterial::matMtrDictionaryFileName() const { return findFile(materialsBaseName_); }

}
