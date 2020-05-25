#ifndef CaloPulseShape_HH
#define CaloPulseShape_HH

// Calculate the values of the digitized pulse shape as a function of the hit time.
// The algorithm is as follow:
//
//   - Import the originl pulse shape, normalize and rebin with the required precision. 
//   - Determine the "start" time of the pulse shape (t0), defined as position where 
//     amplitude = peak amplitude / 1000. Similar for "end" of pulse
//   - for each time bin, calculate the integral of the waveform over the digitized bin width 
//   - extract the digitized pulse from the cached integral values given a start time t0.
//  
// Note 1: Shifting the waveform forward by a time t1 is equivalent to shifting the time 
//         origin backward by a time t1 and leaving the waveform untouched. 
//
// Note 2: We need to build shape after the calibration conditions have been initialized
//

#include <vector>

namespace mu2e {

    class CaloPulseShape
    {
        public:
           CaloPulseShape(double digiSampling, int pulseIntegralSteps, bool doIntegral=false);
           ~CaloPulseShape(){};

           void buildShapes();

	   const std::vector<double>& digitizedPulse(double hitTime) const; 
	   void                       diag() const;

       private:      
	  double                      digiSampling_;
	  int                         pulseIntegralSteps_;	 
	  bool                        doIntegral_;	 
	  int                         nBinShape_;
	  std::vector<double>         integralVal_;
	  mutable std::vector<double> digitizedPulse_;

    };
}
#endif
