// Ed Callaghan
// Interface for a yes/no decision about an individual SimParticle
// November 2024

#ifndef EventMixing_SimParticleSelectionTool_hh
#define EventMixing_SimParticleSelectionTool_hh

#include "Offline/MCDataProducts/inc/SimParticle.hh"

namespace mu2e{
  class SimParticleSelectionTool{
    public:
      SimParticleSelectionTool() = default;
     ~SimParticleSelectionTool() = default;

      virtual bool Select(const SimParticle&) = 0;

    protected:
      /**/

    private:
      /**/
  };
} // namespace mu2e

#endif
