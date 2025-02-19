// Ed Callaghan
// A summable of multiple signals
// February 2025

#include "Offline/TrackerMC/inc/SummedWireSignal.hh"

namespace mu2e{
  SummedWireSignal::SummedWireSignal(): AnalogWireSignal(0.0, 0.0){
    /**/
  }

  SummedWireSignal SummedWireSignal::operator+ (const AnalogWireSignalPtr& rhs){
    this->_components.push_back(rhs);
    return *this;
  }

  double SummedWireSignal::evaluate_shape(double t){
    double rv = 0.0;
    for (auto component: _components){
      rv += component->Evaluate(t);
    }
    return rv;
  }
} // namespace mu2e
