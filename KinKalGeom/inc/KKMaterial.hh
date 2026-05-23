#ifndef KinKalGeom_KKMaterial_hh
#define KinKalGeom_KKMaterial_hh
//
//  build KinKal DetMaterial objects from art parameter configuration
//
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Tuple.h"
// KinKal
#include "KinKal/MatEnv/MatDBInfo.hh"
#include "KinKal/MatEnv/FileFinderInterface.hh"
// KKGeom
#include "Offline/KinKalGeom/inc/KKStrawMaterial.hh"
// mu2e
#include "Offline/Mu2eInterfaces/inc/Detector.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"

#include <memory>
#include <string>

namespace mu2e {
  class KKMaterial : public MatEnv::FileFinderInterface, public Detector {
    public:
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      using MatDBInfo = MatEnv::MatDBInfo;
      struct Config {
        fhicl::Atom<std::string> isotopes { Name("isotopes"), Comment("Filename for istotopes information")};
        fhicl::Atom<std::string> elements { Name("elements"), Comment("Filename for elements information") };
        fhicl::Atom<std::string> materials { Name("materials"), Comment("Filename for materials information") };
        fhicl::Atom<std::string> strawGasMaterialName{ Name("strawGasMaterialName"), Comment("strawGasMaterialName") };
        fhicl::Atom<std::string> strawWallMaterialName{ Name("strawWallMaterialName"), Comment("strawWallMaterialName") };
        fhicl::Atom<std::string> strawWireMaterialName{ Name("strawWireMaterialName"), Comment("strawWireMaterialName") };
        fhicl::Atom<std::string> IPAMaterialName{ Name("IPAMaterialName"), Comment("IPA MaterialName") };
        fhicl::Atom<std::string> STMaterialName{ Name("STMaterialName"), Comment("Stopping Target MaterialName") };
        fhicl::Atom<int> elossMode { Name("IonizationEnergyLossMode"), Comment( "Ionization energy loss mode") };
        fhicl::Atom<double> solidScatter{ Name("SolidScatteringFraction"), Comment("DahlLynch Scattering model cutoff Fraction for solids") };
        fhicl::Atom<double> gasScatter{ Name("GasScatteringFraction"), Comment("DahlLynch Scattering model cutoff Fraction for gases") };
        fhicl::Atom<double> eBrehms{ Name("ElectronBrehmsFraction"), Comment("Electron Brehmsstrahlung cutoff Fraction") };
      };

      explicit KKMaterial( Config const& config, Tracker const& tracker);
      KKStrawMaterial const& strawMaterial() const { return *smat_; }
      auto IPAMaterial() const { return matdbinfo_->findDetMaterial(ipamatname_); }
      auto STMaterial() const { return matdbinfo_->findDetMaterial(stmatname_); }

      // FileFinder interface
      std::string matElmDictionaryFileName() const override;
      std::string matIsoDictionaryFileName() const override;
      std::string matMtrDictionaryFileName() const override;
      std::string findFile( std::string const& path ) const override;
    private:
      // material description files base names (not path)
      std::string elementsBaseName_;
      std::string isotopesBaseName_;
      std::string materialsBaseName_;
      // specific material names
      std::string wallmatname_, gasmatname_, wirematname_,ipamatname_, stmatname_;
      std::unique_ptr<MatDBInfo> matdbinfo_; // material database
      std::unique_ptr<KKStrawMaterial> smat_; // straw material
  };
}
#endif
