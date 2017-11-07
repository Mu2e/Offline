
#ifndef GeneralUtilities_DigitalFiltering_hh
#define GeneralUtilities_DigitalFiltering_hh

#include <vector>

namespace mu2e {
  namespace DigitalFiltering {
    void zpk2tf(std::vector<double> &b, std::vector<double> &a, std::vector<double> &za, std::vector<double> &pa);
    unsigned int iter_factorial(unsigned int n);
    double comb(double n, double k);
    void bilinear(std::vector<double> &bprime, std::vector<double> &aprime, std::vector<double> &b, std::vector<double> &a, double fs);
  }
}

#endif
