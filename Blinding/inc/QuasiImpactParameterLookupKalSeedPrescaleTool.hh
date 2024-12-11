// Ed Callaghan
// art tool to prescale KalSeeds via a lookup-table on the ~minimum DS radial coordinate
// September 2024

#ifndef Blinding_QuasiImpactParameterLookupKalSeedPrescaleTool_hh
#define Blinding_QuasiImpactParameterLookupKalSeedPrescaleTool_hh

// stl
#include <cmath>
#include <vector>

// art
#include "art/Utilities/ToolConfigTable.h"
#include "art/Utilities/ToolMacros.h"

// fhiclcpp
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"
#include "fhiclcpp/types/Sequence.h"

// mu2e
#include "Offline/Blinding/inc/UnivariateLookupKalSeedPrescaleTool.hh"
#include "Offline/DataProducts/inc/SurfaceId.hh"

namespace mu2e{
  class QuasiImpactParameterLookupKalSeedPrescaleTool: public UnivariateLookupKalSeedPrescaleTool{
    public:
      struct Config{
        fhicl::Sequence<std::string> surface_ids{
          fhicl::Name("SurfaceIds"),
          fhicl::Comment("Prioritized sequence of mu2e::SurfaceIds at which tracks may be sampled")
        };
        fhicl::Atom<double> radius_start{
          fhicl::Name("radius_start"),
          fhicl::Comment("Leftmost defined radius value of lookup table")
        };
        fhicl::Atom<double> radius_step_size{
          fhicl::Name("radius_step_size"),
          fhicl::Comment("Abscissa step size of lookup table")
        };
        fhicl::Sequence<double> prescale{
          fhicl::Name("prescale"),
          fhicl::Comment("Sequence containing the coordinates of lookup table")
        };
      };

      using Parameters = art::ToolConfigTable<Config>;
      QuasiImpactParameterLookupKalSeedPrescaleTool(const Parameters&);

    protected:
      SurfaceIdCollection _surface_ids;
      double calculate_observable(const KalSeed&);

    private:
      /**/
  };
} // namespace mu2e

#endif
