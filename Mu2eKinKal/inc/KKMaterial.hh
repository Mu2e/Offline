#ifndef Mu2eKinKal_KKMaterial_hh
#define Mu2eKinKal_KKMaterial_hh
//
//  build KinKal fit configuration objects from art parameter configuration
//
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Tuple.h"
// KinKal
#include "KinKal/MatEnv/MatDBInfo.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawMaterial.hh"
#include "Offline/Mu2eKinKal/inc/KKFileFinder.hh"

#include <memory>
#include <string>

namespace mu2e {
  class KKMaterial {
    public:
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      using MatDBInfo = MatEnv::MatDBInfo;
      struct Config {
        fhicl::Atom<std::string> isotopes { Name("isotopes"), Comment("Filename for istotopes information")};
        fhicl::Atom<std::string> elements { Name("elements"), Comment("Filename for elements information") };
        fhicl::Atom<std::string> materials { Name("materials"), Comment("Filename for materials information") };
        fhicl::Atom<int> eloss { Name("ELossMode"), Comment("Energy Loss model (0=MPV, 1=Moyal"),MatEnv::DetMaterial::moyalmean };
        fhicl::Atom<std::string> strawGasMaterialName{ Name("strawGasMaterialName"), Comment("strawGasMaterialName") };
        fhicl::Atom<std::string> strawWallMaterialName{ Name("strawWallMaterialName"), Comment("strawWallMaterialName") };
        fhicl::Atom<std::string> strawWireMaterialName{ Name("strawWireMaterialName"), Comment("strawWireMaterialName") };
        fhicl::Atom<std::string> IPAMaterialName{ Name("IPAMaterialName"), Comment("IPA MaterialName") };
        fhicl::Atom<std::string> STMaterialName{ Name("STMaterialName"), Comment("Stopping Target MaterialName") };
        fhicl::Atom<double> dahlLynchScatteringFraction{ Name("dahlLynchScatteringFraction"), Comment("dahlLynchScatteringFraction") };
      };

      explicit KKMaterial( Config const& config);
      KKStrawMaterial const& strawMaterial() const;
      auto IPAMaterial() const { return matdbinfo_->findDetMaterial(ipamatname_); }
      auto STMaterial() const { return matdbinfo_->findDetMaterial(stmatname_); }
    private:
      KKFileFinder filefinder_; // used to find material info
      std::string wallmatname_, gasmatname_, wirematname_,ipamatname_, stmatname_;
      MatEnv::DetMaterial::energylossmode eloss_;
      mutable std::unique_ptr<MatDBInfo> matdbinfo_; // material database
      mutable std::unique_ptr<KKStrawMaterial> smat_; // straw material
  };
}
#endif
