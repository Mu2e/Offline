// Bundle together references to various StrawDigi products, to circumvent
// triplicated-loops which pop up when iterating over digis
//
// Ed Callaghan, 2023

#ifndef TrackerMC_StrawDigiBundle_hh
#define TrackerMC_StrawDigiBundle_hh

// mu2e
#include <Offline/RecoDataProducts/inc/StrawDigi.hh>
#include <Offline/MCDataProducts/inc/StrawDigiMC.hh>

namespace mu2e{
  class StrawDigiBundle{
    public:
      StrawDigiBundle(const StrawDigi digi, const StrawDigiADCWaveform adcs, const StrawDigiMC mc);
      StrawDigiBundle(const StrawDigi digi, const StrawDigiADCWaveform adcs);
      StrawDigiBundle(const StrawDigiBundle& bundle);

      const StrawDigi& GetStrawDigi() const;
      const StrawDigiADCWaveform& GetStrawDigiADCWaveform() const;
      const StrawDigiMC& GetStrawDigiMC() const;

      // interface for sorting into buckets of overlapping digitization windows
      const double time() const;
    protected:
      StrawDigi _digi;
      StrawDigiADCWaveform _adcs;
      StrawDigiMC _mc;
    private:
      /**/
  };
}

#endif
