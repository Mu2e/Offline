//
// A spectrum peaked at a given value
//
// C++ includes
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <memory>

// CLHEP includes
#include "CLHEP/Random/RandFlat.h"
#include "cetlib/pow.h"

//ROOT includes
#include "Math/SpecFunc.h"

//Mu2e includes
#include "Mu2eUtilities/inc/PeakedSpectrum.hh"

using cet::pow;
using cet::square;

namespace mu2e {

  PeakedSpectrum::PeakedSpectrum(const fhicl::ParameterSet& psphys, CLHEP::RandFlat& rnFlat):
    _rnFlat (rnFlat                                      ),
    _verbose(psphys.get<int>("peakedSpectrum.verbose", 0))
  {
    _xmax = psphys.get<double>("peakedSpectrum.xmax");
    _xmin = psphys.get<double>("peakedSpectrum.xmin");
    if(_xmin >= _xmax) {
        throw cet::exception("BADCONFIG")
          << "PeakedSpectrum::" <<  __func__ << ": X_min (" <<_xmin << ") >= X_max (" <<_xmax << ")\n";
    }

    std::string spectrum = psphys.get<std::string>("peakedSpectrum.spectrum");
    if(spectrum == "PowerLaw") {
      _spectrum = kPowerLaw;
      //PDF ~ 1 / ((x-a)^b + c)
      _paramsD[0] = psphys.get<double>("peakedSpectrum.peak");   //a
      _paramsI[0] = psphys.get<int>("peakedSpectrum.power");     //b
      _paramsD[1] = psphys.get<double>("peakedSpectrum.offset"); //c
      //normalize the PDF in the given range
      double val1, val2;
      double xrel1 = (_xmax - _paramsD[0]);
      double xrel2 = (_xmin - _paramsD[0]);
      if(_paramsD[1] == 0.) {
        if(_paramsI[0] == 1) {
          val1 = log(xrel1/xrel2); //define difference as ratio inside the log
          val2 = 0.;
        } else {
          val1 = pow(xrel1, 1-_paramsI[0])/(1-_paramsI[0]);
          val2 = pow(xrel2, 1-_paramsI[0])/(1-_paramsI[0]);
        }
      } else {
        //integral in terms of Gauss' hypergeometric function
        val1 = xrel1/_paramsD[1]*ROOT::Math::hyperg(1,1./_paramsI[0],1.+1./_paramsI[0],-1.*pow(xrel1,_paramsI[0])/_paramsD[1]);
        val2 = xrel2/_paramsD[1]*ROOT::Math::hyperg(1,1./_paramsI[0],1.+1./_paramsI[0],-1.*pow(xrel2,_paramsI[0])/_paramsD[1]);
      }
      _norm = 1./(val1-val2);
    } else if(spectrum == "TriangleFlat") {
      _spectrum = kTriangleFlat;
      //PDF ~ max(c*(1-abs(x-a)/b)+d, d)
      _paramsD[0] = psphys.get<double>("peakedSpectrum.peak");          //a
      _paramsD[1] = psphys.get<double>("peakedSpectrum.baseHalfWidth"); //b
      _paramsD[2] = psphys.get<double>("peakedSpectrum.peakHeight");    //c
      _paramsD[3] = psphys.get<double>("peakedSpectrum.flatHeight");    //d
      double xrel1 = (_xmax - _paramsD[0]);
      double xrel2 = (_xmin - _paramsD[0]);
      double integral = _paramsD[3]*(_xmax-_xmin); //flat area
      double triangle_half_area = 0.5*_paramsD[1]*_paramsD[2];
      double xstart = (xrel2 < -_paramsD[1]) ? -_paramsD[1] : xrel2;
      double xend   = (xrel1 >  _paramsD[1]) ?  _paramsD[1] : xrel1;
      if(xstart < 0.) {
        integral += triangle_half_area - 0.5*(_paramsD[1]+xstart); //add left triangle area minus region not included (if any)
        if(xend < 0.) {integral -= triangle_half_area - 0.5*(_paramsD[1]+xend);} //subtrack off end region if not included
      }
      if(xend > 0.) {
        integral += 0.5*(_paramsD[1]-xend); //add right triangle area of included region
        if(xstart > 0.) {integral -= 0.5*(_paramsD[1]-xstart);} //subtrack off start region if not included
      }
      _norm = 1./integral;
    } else {
        throw cet::exception("BADCONFIG")
          << "PeakedSpectrum::" << __func__ << ": Unknown spectrum option " << spectrum.c_str() << "\n";
    }
    if(_verbose > 0) {
      std::cout << "PeakedSpectrum::" << __func__ << ": " << _xmin << " < x < " << _xmax
                << ", using Spectrum " << spectrum.c_str() << " with norm = " << _norm << std::endl;
    }
  }


  double PeakedSpectrum::getWeight(double x) const {
    double weight(0.);
    if(_spectrum == kPowerLaw) {
      weight = 1./(pow(x-_paramsD[0],_paramsI[0]) + _paramsD[1]);
    } else if(_spectrum == kTriangleFlat) {
      weight = std::max(_paramsD[2]*(1.-std::abs(x-_paramsD[0])/_paramsD[1]), _paramsD[3]);
    }
    weight *= _norm;
    if(_verbose > 9) {
      std::cout << "PeakedSpectrum::" << __func__ << ": x = " << x << " weight = " << weight << std::endl;
    }
    return weight;
  }

  double PeakedSpectrum::getMax() const {
    return _paramsD[0]; //peak value in all spectra
  }

//-----------------------------------------------------------------------------
  double PeakedSpectrum::fire() const {
    double pdfMax = getMax();
    double x(0.), prob(0.), threshold(0.);
    do {
      x         = _xmin + (_xmax-_xmin)*_rnFlat.fire();
      threshold = getWeight(x);
      prob      = pdfMax*_rnFlat.fire();
    } while (prob > threshold);
    return x;
  }
}
