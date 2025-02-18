// Ed Callaghan
// Interface to define an analog signal, which is just a function of time
// February 2025

#ifndef TrackerMC_AnalogWireSignal_hh
#define TrackerMC_AnalogWireSignal_hh

// stl
#include <memory>

// boost
#include "boost/math/tools/roots.hpp"

// mu2e
#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/DataProducts/inc/TrkTypes.hh"
#include "Offline/TrackerConditions/inc/StrawElectronics.hh"

namespace mu2e{
  class AnalogWireSignal{
    public:
      AnalogWireSignal(double, double);
      ~AnalogWireSignal() = default;

      virtual double Evaluate(double) = 0;
      virtual bool CrossesThreshold(double, double);
      virtual bool CoarseThresholdCrossingTime(double, double, double&);
      virtual double ThresholdCrossingTime(double, double, double, double);

      void DigitalTimeOverThreshold(const StrawElectronics&,
                                      const double,
                                      const double,
                                      TrkTypes::TOTValue&);
      void Digitize(const StrawElectronics&,
                    const StrawId&,
                    const double,
                    TrkTypes::ADCTimes&,
                    TrkTypes::ADCWaveform&,
                    TrkTypes::ADCValue&);

    protected:
      double _time_lo; // lower window bound
      double _time_hi; // upper window bound
    private:
      /**/
  };

  typedef std::shared_ptr<AnalogWireSignal> AnalogWireSignalPtr;
} // namespace mu2e

#endif
