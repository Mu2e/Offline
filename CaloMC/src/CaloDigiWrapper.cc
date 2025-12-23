// Ed Callaghan
// Simple wrapper around calorimeter digis to interface into a partitioning scheme
// September 2025

#include "Offline/CaloMC/inc/CaloDigiWrapper.hh"

namespace mu2e{
  CaloDigiWrapper::CaloDigiWrapper(const CaloDigi& digi):
      _digi(digi.SiPMID(), digi.t0(), digi.waveform(), digi.peakpos()){
    /**/
  }

  const CaloDigi& CaloDigiWrapper::Digi() const{
    const auto& rv = _digi;
    return rv;
  }

  const double CaloDigiWrapper::time() const{
    auto t0 = _digi.t0();
    auto rv = static_cast<double>(t0);
    return rv;
  }
} // namespace mu2e
