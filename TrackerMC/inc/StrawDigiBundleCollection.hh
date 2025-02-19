// Containerize StrawDigiBundles, with a extra functionality for resolving
// situations where multiple digis have overlapping digitization windows
//
// Ed Callaghan, 2023

#ifndef TrackerMC_StrawDigiBundleCollection_hh
#define TrackerMC_StrawDigiBundleCollection_hh

// stl
#include <map>
#include <memory>
#include <vector>

// art
#include "art/Framework/Principal/Handle.h"

// cetlib_except
#include "cetlib_except/exception.h"

// mu2e
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/TrackerConditions/inc/StrawElectronics.hh"
#include "Offline/TrackerMC/inc/StrawDigiBundle.hh"

namespace mu2e{
  using SDBC_iterator = std::vector<StrawDigiBundle>::iterator;
  using SDBC_const_iterator = std::vector<StrawDigiBundle>::const_iterator;

  class StrawDigiBundleCollection{
    public:
      // forwarded calls to underlying container
      size_t size() const;
      SDBC_iterator begin();
      SDBC_const_iterator begin() const;
      SDBC_iterator end();
      SDBC_const_iterator end() const;
      StrawDigiBundle& operator[](size_t i);
      const StrawDigiBundle& operator[](size_t i) const;

      // convenience methods for accepting new StrawDigiBundles
      // from different source situations
      void Append(const StrawDigiCollection& digis,
                  const StrawDigiADCWaveformCollection& adcss,
                  const StrawDigiMCCollection& mcs);
      void Append(const StrawDigiCollection& digis,
                  const StrawDigiADCWaveformCollection& adcss);
      void Append(const StrawDigiBundle bundle);

      // utility functions
      void FillStrawDigis(StrawDigiCollection& rv);
      void FillStrawDigiADCWaveforms(StrawDigiADCWaveformCollection& rv);
      void FillStrawDigiMCs(StrawDigiMCCollection& rv);
      StrawDigiBundleCollection operator+= (const StrawDigiBundleCollection& other);

      // accessors, as objects
      StrawDigiCollection GetStrawDigis();
      StrawDigiADCWaveformCollection GetStrawDigiADCWaveforms();
      StrawDigiMCCollection GetStrawDigiMCs();

      // accessors, as pointers
      std::unique_ptr<StrawDigiCollection> GetStrawDigiPtrs();
      std::unique_ptr<StrawDigiADCWaveformCollection> GetStrawDigiADCWaveformPtrs();
      std::unique_ptr<StrawDigiMCCollection> GetStrawDigiMCPtrs();

      // identify sets of  digis with overlapping digitization windows, and
      // reduce each such set to a single digi, representing their "sum"
      void ResolveCollisions(const StrawElectronics& conditions, StrawDigiBundleCollection& rv);
    protected:
      std::vector<StrawDigiBundle> _bundles;
      void ResolveCollision(StrawDigiBundleCollection& collided, StrawDigiBundleCollection& rv);
    private:
      /**/
  };
}

#endif
