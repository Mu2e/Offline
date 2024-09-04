// Ed Callaghan
// art tool to prescale KalSeeds via a lookup-table on the momentum
// September 2024

#ifndef Blinding_MomentumLookupKalSeedPrescaleTool_hh
#define Blinding_MomentumLookupKalSeedPrescaleTool_hh

// stl
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
  class MomentumLookupKalSeedPrescaleTool: public UnivariateLookupKalSeedPrescaleTool{
    public:
      struct Config{
        fhicl::Sequence<std::string> surface_ids{
          fhicl::Name("SurfaceIds"),
          fhicl::Comment("Prioritized sequence of mu2e::SurfaceIds at which tracks may be sampled")
        };
        fhicl::Atom<double> momentum_start{
          fhicl::Name("momentum_start"),
          fhicl::Comment("Leftmost defined momentum value of lookup table")
        };
        fhicl::Atom<double> momentum_step_size{
          fhicl::Name("momentum_step_size"),
          fhicl::Comment("Abscissa step size of lookup table")
        };
        fhicl::Sequence<double> prescale{
          fhicl::Name("prescale"),
          fhicl::Comment("Sequence containing the coordinates of lookup table")
        };
      };

      using Parameters = art::ToolConfigTable<Config>;
      MomentumLookupKalSeedPrescaleTool(const Parameters&);

    protected:
      SurfaceIdCollection _surface_ids;
      double calculate_observable(const KalSeed&);

    private:
      /**/
  };
} // namespace mu2e

#endif
