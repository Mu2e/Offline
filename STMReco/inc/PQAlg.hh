#ifndef PQAlg_hh
#define PQAlg_hh

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <math.h>
#include <time.h>
#include "cetlib_except/exception.h"

#include "Offline/RecoDataProducts/inc/STMDigi.hh"

namespace mu2e {

  class PQAlg {

  public:

    //Constructors
    PQAlg() {
      //All the parameters value are just temporary, needs to be optimized. Currently set at 4ns/pt
      delta_PF=16;
      delta_Walkback=8;
      delta_Integration=64;
      Gap=8;
      delta_PB=16;
      total_length=delta_PB+delta_PF+delta_Integration+Gap;
      N_shift=5;
      a0=9,b0=10;
      i=0; j=0; //numbers for counting iterations on loops; see further explanantion when being used
      bound1=25; bound2=50; //numerical bounds for assesing the fron_back pedestal difference affecting pulse integration; see further explanantion when being used
    }

    //PQAlg methods

    /* The pulse energy calculating function, which reads from a big txt file, each line contains a floating type data as the instantaneous value of the pulse voltage. each pulse is in a 'pulse segment'. The whole file contains N_pulse 'pulse segments'. Currently, each pulse segment has 502 floating points, and there are 100x700=70,000 pulse segments in total.

       parameters explained; All the parameters being used are based on the FADC voltage data points are sampled at 500MHz(2ns/4ns) interval

       Inputs
       buffer_length:           same as total_length
       Baseline_level_bounds:   check the  difference in pulse integration when subtracting the back_pedestal rather than the                          front_pedestal.
       bound1:25(the difference=4% of pulse integration); bound2:50(the difference=2% of pulse integration)
       Outputs

       pulse_time:             a reference point of each pulse to indicates its arriving time, not used in energy calculation
       test_y:                  store pulse integration value of each pulse, proportional to pulse energy
       Quality factor:          a discrete value outputs 0,1 or 2 describing the pulse integration quality

    */
    int process_pulse(const mu2e::STMDigi& digi);

    size_t getNFound() { return _energies.size(); }
    const int& getQuality(size_t i) const { return _qualities.at(i); }
    const double& getEnergy(size_t i) const { return _energies.at(i); }
    const int& getTrigSample(size_t i) const { return _trigSamples.at(i); }

  private:
    // parameters will go here
    int delta_PF; // length of time of front pedestal(baseline), expressed in number of data points taken
    int delta_Walkback; // back off time after a trigger is found, before integrating the pulse
    int delta_Integration; // Pulse integration time, starting from the walkback point;
    int Gap; // time between integration period and back pedestal
    int delta_PB; // length of time of back pedestal(baseline)
    int total_length; // data sample length in one pulse
    int N_shift; // shift reverse pulse by N_shift intervals as to the normal pulse to perform the CFD pulse time method
    int a0,b0; // relative height of the normal and reversed pulse in CFD timing method

    int i; int j; //numbers for counting iterations on loops; see further explanantion when being used
    int bound1; int bound2; //numerical bounds for assesing the fron_back pedestal difference affecting pulse integration; see further explanantion when being used

    std::vector<double> _energies;
    std::vector<int> _qualities;
    std::vector<int> _trigSamples;
  };
};
#endif
