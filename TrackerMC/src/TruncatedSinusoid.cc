// Ed Callaghan
// A fixed-length sinusoid
// February 2025

// mu2e
#include "Offline/TrackerMC/inc/TruncatedSinusoid.hh"

namespace mu2e{
  TruncatedSinusoid::TruncatedSinusoid(double amplitude,
                                       double frequency,
                                       double phase,
                                       double time_lo,
                                       double time_hi):
      _amplitude(amplitude),
      _frequency(frequency),
      _phase(phase),
      _time_lo(time_lo),
      _time_hi(time_hi){
    /**/
  }

  double TruncatedSinusoid::Evaluate(double t){
    double rv = 0.0;
    if ((_time_lo < t) && (t < _time_hi)){
      double phase = CLHEP::twopi * t * _frequency;
      rv = _amplitude * sin(phase - _phase);
    }
    return rv;
  }
} // namespace mu2e
