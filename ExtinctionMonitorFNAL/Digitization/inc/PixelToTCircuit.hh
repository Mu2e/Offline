// Simulate pulse shape in the pixel time over threshold (ToT) circuit.
//
// Andrei Gaponenko, 2012

#ifndef ExtinctionMonitorFNAL_Digitization_PixelToTCircuit_hh
#define ExtinctionMonitorFNAL_Digitization_PixelToTCircuit_hh

namespace mu2e {
  namespace ExtMonFNAL {


    class PixelToTCircuit {
    public:

      // For the given threshold  and qCalib calibration injected charge
      // the calibrated circuit should produce average ToT = totCalib
      PixelToTCircuit(double threshold, double qCalib, int totCalib, double clockTick);

      void addCharge(double q) { q_ += q; }

      // get current charge
      double charge() const { return q_; }
      bool   high() const { return q_ > threshold_; }

      // compute discharge time to go below threshold
      double computeTrailingEdge() const;

      // discharge the capasitor
      void wait(double time);

    private:
      double q_; // current charge
      double threshold_;
      double k_; // the discharge slope
    };

  } // namespace ExtMonFNAL
} // namespace mu2e

#endif
