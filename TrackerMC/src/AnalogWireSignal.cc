// Ed Callaghan
// Interface to define an analog signal, which is just a function of time
// February 2025

// mu2e
#include "Offline/TrackerMC/inc/AnalogWireSignal.hh"

namespace mu2e{
  AnalogWireSignal::AnalogWireSignal(UnaryFunctionPtr shape): _shape(shape){
  }

  // evaluate the signal by shifting the primary shape, and adding the rest
  double AnalogWireSignal::Evaluate(double time){
    double rv = _shape->Evaluate(time - _delay);
    for (auto summand: _summands){
      rv += summand.Evaluate(time);
    }
    return rv;
  }

  void AnalogWireSignal::AddDelay(double delay){
    _delay += delay;
  }

  AnalogWireSignal AnalogWireSignal::operator+ (const AnalogWireSignal& rhs){
    _summands.push_back(rhs);
    return *this;
  }

  // convenience wrapper
  bool AnalogWireSignal::CrossesThreshold(double threshold,
                                          double time_lo,
                                          double time_hi,
                                          double spacing){
    double tmp;
    bool rv = this->CoarseThresholdCrossingTime(threshold, time_lo, time_hi, spacing, tmp);
    return rv;
  }

  // default coarse approach to find threshold-crossing via linear scan
  bool AnalogWireSignal::CoarseThresholdCrossingTime(double threshold,
                                                     double time_lo,
                                                     double time_hi,
                                                     double spacing,
                                                     double& out){
    for (double time = time_lo ; time < time_hi ; time += spacing){
      if (threshold < this->Evaluate(time)){
        out = time;
        return true;
      }
    }
    return false;
  }

  // default fine approach to determine threshold-crossing via bisection method
  // the boost-backed implementation brackets the crossing, assuming it exists,
  // to within some gross time margin. the root is approximated as the midpoint
  // of this bracket
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

  // determine if a threshold-crossing exists and, if so, delay the primary
  // signal s.t. the crossing occurs at a predetermined time
  bool AnalogWireSignal::TranslateToThresholdCrossingTime(double threshold,
                                                          double time,
                                                          double time_lo,
                                                          double time_hi,
                                                          double coarse_spacing,
                                                          double tolerance){
    double crTime;
    bool crossed = this->CoarseThresholdCrossingTime(threshold,
                                                     time_lo,
                                                     time_hi,
                                                     coarse_spacing,
                                                     crTime);
    if (crossed){
      crTime = this->ThresholdCrossingTime(threshold,
                                           crTime - coarse_spacing,
                                           crTime + coarse_spacing,
                                           tolerance);
      this->AddDelay(time - crTime);
    }
    return crossed;
  }

  // compute digital TOT value
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

  // compute digital waveform and amplitude around an assumed threshold-crossing
  void AnalogWireSignal::Digitize(const StrawElectronics& electronics,
                                  const StrawId& sid,
                                  const double crTime,
                                  TrkTypes::ADCTimes& sample_times,
                                  TrkTypes::ADCWaveform& sample_values,
                                  TrkTypes::ADCValue& pmp){
    // calculate sample times
    electronics.adcTimes(crTime, sample_times);

    // compute saturated analog signal
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
