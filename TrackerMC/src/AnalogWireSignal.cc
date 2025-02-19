// Ed Callaghan
// Interface to define an analog signal, which is just a function of time
// February 2025

// mu2e
#include "Offline/TrackerMC/inc/AnalogWireSignal.hh"

namespace mu2e{
  AnalogWireSignal::AnalogWireSignal(double time_lo, double time_hi):
      _time_lo(time_lo),
      _time_hi(time_hi){
    /**/
  }

  double AnalogWireSignal::Evaluate(double time){
    double rv = this->evaluate_shape(time);
    return rv;
  }

  bool AnalogWireSignal::CrossesThreshold(double threshold, double spacing){
    double tmp;
    bool rv = this->CoarseThresholdCrossingTime(threshold, spacing, tmp);
    return rv;
  }

  bool AnalogWireSignal::CoarseThresholdCrossingTime(double threshold,
                                                     double spacing,
                                                     double& out){
    for (double time = _time_lo ; time < _time_hi ; time += spacing){
      if (threshold < this->Evaluate(time)){
        out = time;
        return true;
      }
    }
    return false;
  }

  double AnalogWireSignal::ThresholdCrossingTime(double threshold,
                                                 double lhs,
                                                 double rhs,
                                                 double tolerance){
    auto f = [this, threshold] (double time) {
      double rv = this->Evaluate(time) - threshold;
      return rv;
    };
    auto tol = [tolerance] (double x, double y){
      bool rv = fabs(x - y) < tolerance;
      return rv;
    };
    const auto pair = boost::math::tools::bisect(f, lhs, rhs, tol);
    double rv = 0.5 * (pair.first + pair.second);
    return rv;
  }

  void AnalogWireSignal::DigitalTimeOverThreshold(
                                          const StrawElectronics& electronics,
                                          const double threshold,
                                          const double crTime,
                                          TrkTypes::TOTValue& tot){
    auto lsb = electronics.totLSB();
    auto max = electronics.maxTOT();
    tot = max;

    double effective_threshold = threshold - electronics.triggerHysteresis();
    for (size_t i = 0 ; i < max ; i++){
      double time = crTime + i*lsb;
      double value = this->Evaluate(time);
      if (value < effective_threshold){
        tot = i;
        return;
      }
    }
  }

  void AnalogWireSignal::Digitize(const StrawElectronics& electronics,
                                  const StrawId& sid,
                                  const double crTime,
                                  TrkTypes::ADCTimes& sample_times,
                                  TrkTypes::ADCWaveform& sample_values,
                                  TrkTypes::ADCValue& pmp){
    // calculate sample times
    electronics.adcTimes(crTime, sample_times);

    // calculate saturated analog signal
    TrkTypes::ADCVoltages voltages(sample_times.size());
    for (size_t i = 0 ; i < sample_times.size() ; i++){
      double voltage = this->Evaluate(sample_times[i]);
      double saturated = electronics.saturatedResponse(voltage);
      voltages[i] = saturated;
    }

    // digitize
    electronics.digitizeWaveform(sid, voltages, sample_values, pmp);
  }
} // namespace mu2e
