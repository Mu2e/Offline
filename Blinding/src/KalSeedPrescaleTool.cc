// Ed Callaghan
// Interface for art tool to calculate a KalSeed-dependent functional prescale
// September 2024

#include "Offline/Blinding/inc/KalSeedPrescaleTool.hh"

namespace mu2e{
  double KalSeedPrescaleTool::AcceptanceRate(const KalSeed& kalseed){
    auto rv = this->calculate_acceptance_rate(kalseed);
    return rv;
  }
} // namespace mu2e
