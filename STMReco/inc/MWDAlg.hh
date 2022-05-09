#ifndef MWDAlg_hh
#define MWDAlg_hh

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <utility> // std::pair
#include <stdexcept> // std::runtime_error
#include <sstream> // std::stringstream
#include <fstream>
#include <boost/chrono.hpp>
#include <cmath>
#include <sys/stat.h>
#include <cstring>

namespace mu2e {

  // some classes that are used in the algorithm
  // TODO: clean up
  class data {

  public:
    const int16_t* adc;
    uint64_t t0;
    uint32_t nadc;
  };

  class stats {
  public :

    static double mean(double* x, int n) {
      double sum = 0.0;
      double N = (double) n;
      for (int i = 0; i < n; i++){
        sum = sum + x[i];
      }
      return sum/(double) N;
    }

    static double mean(std::vector<double> x){
      double sum = 0.0;
      for (unsigned int i = 0; i < x.size(); i++){
        sum = sum + x[i];
      }
      return sum/(double) x.size();
    }

    static double rms(double* x, int n) {
      double mean_val = mean(x,n);
      return rms(x, mean_val, n);
    }

    static double rms(double* x, double mean_val, int n) {
      if (n <= 1) {
        std::cout << "ERROR: stats::rms - trying to calculate rms of array with one or fewer elements : " << n << std::endl;
        exit(-1);
      }
      double variance = 0.0;
      for(int i = 0; i < n; i++){
        variance += pow(x[i] - mean_val, 2);
      }
      double N = (double) n;
      return sqrt(variance/N);
    }

    static double rms(std::vector<double> x) {
      double mean_val = mean(x);
      return rms(x, mean_val);
    }

    static double rms(std::vector<double> x, double mean_val) {
      if (x.size() <= 1) {
        std::cout << "ERROR: stats::rms - trying to calculate rms of array with one or fewer elements : " << x.size() << std::endl;
        exit(-1);
      }
      double variance = 0.0;
      for(unsigned int i = 0; i < x.size(); i++){
        variance += pow(x[i] - mean_val, 2);
      }
      double N = (double) x.size();
      return sqrt(variance/N);
    }


  };

  class peaks {

  public:
    int npeaks;
    uint32_t trigger_number;
    uint16_t trigger_type;
    uint32_t nadc_values;
    std::vector<double> peak_heights;
    std::vector<double> peak_times;
  };



  class MWDAlg {

  public:

    //Constructors
    MWDAlg();
    MWDAlg(double _M, double _L, double _tau, double _nsigma_cut, double _thresholdgrad, double _fADC, int _cut_mode, double _fixed_cut);

    //MWDAlg methods
    std::string print();
    void mwd_algorithm(data* adc_values);
    std::vector<double> calculate_baseline();
    peaks* find_peaks(double baseline_mean, double baseline_rms, double time_offset); // time offset is in us.

  private:
    int M, L, cut_mode;
    double tau, nsigma_cut, thresholdgrad, fADC, fixed_cut_parameter;
    double threshold_cut;

    double* l; // output of MWDAlg algorithm
    int nadc;
  };
};
#endif
