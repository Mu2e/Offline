// Ed Callaghan
// Simple wrapper around calorimeter digis to interface into a partitioning scheme
// September 2025

#ifndef CaloMC_CaloDigiWrapper_hh
#define CaloMC_CaloDigiWrapper_hh

// mu2e
#include "Offline/RecoDataProducts/inc/CaloDigi.hh"

namespace mu2e{
  class CaloDigiWrapper{
    public:
      using SiPMID_t = int;
      using sample_t = int;
      using pos_t = size_t;
      CaloDigiWrapper(const CaloDigi&);
      const CaloDigi& Digi() const;

      // interface for sorting into buckets of overlapping digitization windows
      const double time() const;

    protected:
      const CaloDigi _digi;

    private:
      /**/
  };
} // namespace mu2e

#endif
