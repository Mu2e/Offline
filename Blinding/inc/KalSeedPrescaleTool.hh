// Ed Callaghan
// Interface for art tool to calculate a KalSeed-dependent functional prescale
// September 2024

#ifndef Blinding_KalSeedPrescaleTool_hh
#define Blinding_KalSeedPrescaleTool_hh

#include "Offline/RecoDataProducts/inc/KalSeed.hh"

namespace mu2e{
  class KalSeedPrescaleTool{
    public:
      KalSeedPrescaleTool() = default;
      double AcceptanceRate(const KalSeed&);

    protected:
      virtual double calculate_acceptance_rate(const KalSeed&) = 0;

    private:
      /**/
  };
} // namespace mu2e

#endif
