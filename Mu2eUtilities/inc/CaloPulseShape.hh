#ifndef CaloPulseShape_HH
#define CaloPulseShape_HH

// Calculate the values of the digitized pulse shape as a function of the hit time.
// The value stored are the integral of the waveform over the digitization bin width.
//
// The waveform "starting point" corresponds to the time the PE hit the readout.
// The definition is arbitrary, but it must be internally consistent!
//
// Shifting the waveform forward by a time dt is equivalent to shifting the time origin backward by dt
//  - value at nbin-i correspond to waveform shifted by time nbin+i,
//  - t0 value is located at bin nSteps_
//  - the t0 value correspond to the content of the digitized bin whose LOW EDGE is at time t0
//
// 1) digitizedPulse(hitTime) returns a waveform with hitTime corresponding to low edge of first bin
// 2) evaluate(deltaTime) return value of digitized bin at a given time difference with peak time value
//
//  NOTE: uncomment the pline creation if the discontinuities in the second order derivative arising from the
//        linear piecewise approxmiation are problematic for the minimization

#include <vector>
#include <string>

namespace mu2e {

    class CaloPulseShape
    {
       public:
          CaloPulseShape(const std::string& fileName, const std::string& histName, double digiSampling);
          ~CaloPulseShape() {};

          void buildShapes();

          const std::vector<double>& digitizedPulse  (double hitTime)        const;
          double                     evaluate        (double timeDifference) const;
          double                     fromPeakToT0    (double timePeak)       const;
          void                       diag            (bool fullDiag=false)   const;

       private:
          std::string                 fileName_;
          std::string                 histName_;
          int                         nSteps_;
          double                      digiStep_;
          int                         nBinShape_;
          std::vector<double>         pulseVec_;
          double                      deltaT_;
          mutable std::vector<double> digitizedPulse_;
    };

}
#endif
