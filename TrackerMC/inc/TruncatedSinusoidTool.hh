// Ed Callaghan
// Produce sinusoidal waves
// February 2025

#ifndef TrackerMC_TruncatedSinusoidTool_hh
#define TrackerMC_TruncatedSinusoidTool_hh

#include "Offline/TrackerMC/inc/AnalogSignalShapeTool.hh"
#include "Offline/TrackerMC/inc/TruncatedSinusoid.hh"

namespace mu2e{
  class TruncatedSinusoidTool: public AnalogSignalShapeTool{
    public:
      // note that we do not choose the phase
      // if this wave triggers, then the phase is fixed by the crossing-time
      struct Config{
        fhicl::Atom<double> amplitude{
          fhicl::Name("amplitude"),
          fhicl::Comment("Wave amplitude in milliVolts")
        };
        fhicl::Atom<double> frequency{
          fhicl::Name("frequency"),
          fhicl::Comment("Wave frequency in gigahertz")
        };
        fhicl::Atom<double> phase{
          fhicl::Name("phase"),
          fhicl::Comment("Wave phase in nanoseconds")
        };
        fhicl::Atom<double> time_lo{
          fhicl::Name("time_lo"),
          fhicl::Comment("Start time of signal in nanoseconds")
        };
        fhicl::Atom<double> time_hi{
          fhicl::Name("time_hi"),
          fhicl::Comment("End time of signal in nanoseconds")
        };
      };

      using Parameters = art::ToolConfigTable<Config>;
      TruncatedSinusoidTool(const Parameters&);
     ~TruncatedSinusoidTool() = default;

      virtual UnaryFunctionPtr Sample() override;

    protected:
      double _amplitude;
      double _frequency;
      double _phase;
      double _time_lo;
      double _time_hi;

    private:
      /**/
  };
} // namespace mu2e

#endif
