// Ed Callaghan
// Stop-gap tool to fake a volume lookup, only valid for protonabs1 (IPA)
// November 2024

#ifndef EventMixing_InnerProtonAbsorberPseudoVolumeLookupTool_hh
#define EventMixing_InnerProtonAbsorberPseudoVolumeLookupTool_hh

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
#include "Offline/MECOStyleProtonAbsorberGeom/inc/MECOStyleProtonAbsorber.hh"

namespace mu2e{
  class InnerProtonAbsorberPseudoVolumeLookupTool{
    public:
      struct Config{
        fhicl::Atom<std::string> ipa_name{
          fhicl::Name("IPA"),
          fhicl::Comment("Volume pseudoname to assign to ipa coordinates")
        };
        fhicl::Atom<std::string> other_name{
          fhicl::Name("Other"),
          fhicl::Comment("Volume pseudoname to assign to non-ipa coordinates")
        };
      };

      using Parameters = art::ToolConfigTable<Config>;
      InnerProtonAbsorberPseudoVolumeLookupTool(const Parameters&);
      virtual ~InnerProtonAbsorberPseudoVolumeLookupTool();

      // volume pseudo-lookup
      virtual std::string Volume(const CLHEP::Hep3Vector&);
      virtual std::string StartVolume(const SimParticle&);

    protected:
      std::string _ipa_name;
      std::string _other_name;
      std::unique_ptr< GeomHandle<DetectorSystem> > _frame;
      std::unique_ptr< GeomHandle<MECOStyleProtonAbsorber> > _ipa;

      // deferred initialization is necessary to make use of GeometryService
      bool _initialized;
      void initialize();

    private:
      /**/
  };
} // namespace mu2e

#endif
