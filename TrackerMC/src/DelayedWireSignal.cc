// Ed Callaghan
// A delay added to a preexisting signal
// February 2025

#include "Offline/TrackerMC/inc/DelayedWireSignal.hh"

namespace mu2e{
  DelayedWireSignal::DelayedWireSignal(AnalogWireSignalPtr& signal,
                                       double delay):
      AnalogWireSignal(0.0, 0.0),
      _signal(signal),
      _delay(delay){
    /**/
  }

  double DelayedWireSignal::evaluate_shape(double t){
    t -= _delay;
    double rv = this->_signal->Evaluate(t);
    return rv;
  }
} // namespace mu2e
