// Ed Callaghan
// Interface to define an analog signal, which is just a function of time
// February 2025

#ifndef TrackerMC_AnalogWireSignal_hh
#define TrackerMC_AnalogWireSignal_hh

// stl
#include <vector>

// boost
#include "boost/math/tools/roots.hpp"

// mu2e
#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/DataProducts/inc/TrkTypes.hh"
#include "Offline/GeneralUtilities/inc/UnaryFunction.hh"
#include "Offline/TrackerConditions/inc/StrawElectronics.hh"

namespace mu2e{
  using UnaryFunctionPtr = std::shared_ptr<UnaryFunction>;

  class AnalogWireSignal{
    public:
      AnalogWireSignal(UnaryFunctionPtr);
      ~AnalogWireSignal() = default;

      double Evaluate(double);
      void AddDelay(double);
      AnalogWireSignal operator+ (const AnalogWireSignal&);

      virtual bool CrossesThreshold(double, double, double, double);
      virtual bool CoarseThresholdCrossingTime(double, double, double,
                                               double, double&);
      virtual double ThresholdCrossingTime(double, double, double, double);

      bool TranslateToThresholdCrossingTime(double, double,
                                            double, double,
                                            double, double);

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
      double _delay;
      UnaryFunctionPtr _shape;
      std::vector<AnalogWireSignal> _summands;

    private:
      /**/
  };
} // namespace mu2e

#endif
