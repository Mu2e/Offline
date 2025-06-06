// Ed Callaghan
// Stop-gap tool to fake a volume lookup, only valid for protonabs1 (IPA)
// November 2024

#ifndef EventMixing_PseudoCylindricalVolumeLookupTool_hh
#define EventMixing_PseudoCylindricalVolumeLookupTool_hh

// stl
#include <string>

// art
#include "art/Utilities/ToolConfigTable.h"
#include "art/Utilities/ToolMacros.h"

// clhep
#include "CLHEP/Vector/ThreeVector.h"

// fhiclcpp
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"

// mu2e
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/BeamlineGeom/inc/ProtonAbsorber.hh"
#include "Offline/StoppingTargetGeom/inc/StoppingTarget.hh"
#include "Offline/StoppingTargetGeom/inc/TargetFoil.hh"

namespace mu2e{
  class PseudoCylindricalVolumeLookupTool{
    public:
      struct Config{
        fhicl::Atom<std::string> ipa_name{
          fhicl::Name("IPA"),
          fhicl::Comment("Volume pseudoname to assign to ipa coordinates")
        };
        fhicl::Atom<std::string> st_name{
          fhicl::Name("ST"),
          fhicl::Comment("Volume pseudoname to assign to s.t. coordinates")
        };
        fhicl::Atom<std::string> other_name{
          fhicl::Name("Other"),
          fhicl::Comment("Volume pseudoname to assign to non-ipa coordinates")
        };
      };

      using Parameters = art::ToolConfigTable<Config>;
      PseudoCylindricalVolumeLookupTool(const Parameters&);
      virtual ~PseudoCylindricalVolumeLookupTool();

      // volume pseudo-lookup
      virtual std::string Volume(const CLHEP::Hep3Vector&);
      virtual std::string StartVolume(const SimParticle&);

    protected:
      std::string _ipa_name;
      std::string _st_name;
      std::string _other_name;
      std::unique_ptr< GeomHandle<DetectorSystem> > _frame;
      std::unique_ptr< GeomHandle<ProtonAbsorber> > _ipa;
      std::unique_ptr< GeomHandle<StoppingTarget> > _st;

      // deferred initialization is necessary to make use of GeometryService
      bool _initialized;
      void initialize();

      bool in_cylindrical_shell(const CLHEP::Hep3Vector&,
                                double, double,
                                double, double);

    private:
      /**/
  };
} // namespace mu2e

#endif
