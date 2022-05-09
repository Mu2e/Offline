///////////////////////////////////////
//
// Zero-suppression algorithm for STM
// Original authors: Claudia Alvarez-Garcia, Alex Kesharvarzi, and Mark Lancaster
// Adapated for Offline: Andy Edmonds
//
#ifndef ZPAlg_H
#define ZPAlg_H

#include<iostream>
#include <string>
#include <fstream>
#include <vector>
#include <utility> // std::pair
#include <stdexcept> // std::runtime_error
#include <sstream> // std::stringstream
#include<fstream>
#include <boost/chrono.hpp>
#include<cmath>
#include <sys/stat.h>
#include <cstring>

using namespace std;

namespace mu2e {
  class ZPAlg {
  public:

    //Constructors
    //Default
    ZPAlg(double fadc, double tb, double ta, double thresh) :
      fADC(fadc*1e6), tbefore(tb*1e-3), tafter(ta*1e-3), threshold(thresh){}

    int16_t* ZeroSuppress (const int16_t* ADC, unsigned long int n, std::vector<size_t>& starts, std::vector<size_t>& ends);

  private:

    double fADC; //Hz
    //store 2 microseconds of data to the left of the trigger
    double tbefore;
    //store 20 microseconds of data to the right of the trigger
    double tafter;
    //Find the trigger
    int threshold;
  };
}
#endif
