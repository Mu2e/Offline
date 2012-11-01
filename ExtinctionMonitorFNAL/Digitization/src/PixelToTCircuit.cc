// Simulate pulse shape in the pixel time over threshold (ToT) circuit.
//
// Andrei Gaponenko, 2012

#include "ExtinctionMonitorFNAL/Digitization/inc/PixelToTCircuit.hh"

#include <algorithm>

namespace mu2e {
  namespace ExtMonFNAL {

    PixelToTCircuit::PixelToTCircuit(double threshold, double qCalib, int totCalib, double clockTick)
      : q_(0.), threshold_(threshold), k_()
    {
      k_ = (qCalib - threshold)/(totCalib*clockTick);
    }

    void PixelToTCircuit::wait(double t) {
      q_ = std::max(0., q_ - k_*t);
    }

    double PixelToTCircuit::computeTrailingEdge() const {
      return (q_ - threshold_)/k_;
    }

  } // namespace ExtMonFNAL
} // namespace mu2e
