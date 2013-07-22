#ifndef Mu2eUtilities_BinnedSpectrum_hh
#define Mu2eUtilities_BinnedSpectrum_hh

//
// Base class to allow generic access to all the classes that define
// a momentum spectrum.
//
// $Id: BinnedSpectrum.hh,v 1.1 2013/07/22 18:57:42 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2013/07/22 18:57:42 $
//
// Original author Kyle Knoepfel 
//                 

// C++ includes
#include <iostream>
#include <string>
#include <utility>
#include <vector>

namespace mu2e {

  class BinnedSpectrum {

  public:

    // To make accessible to CLHEP::RandGeneral, return pointer to 
    // first spectrum entry
    const double* getPDF()   const { return &(*_spectrum.begin() ); }
    unsigned      getNbins() const { return _spectrum.size();       }

    template<class Shape, typename... Args>
    void initialize(double emin, double emax, double spectRes, Args... args) {
      
      Shape s( std::forward<Args>(args)... );
      for ( double step = emin ; step <= emax ; step += spectRes ) {
        //        if ( step > emin ) break;
        _spectrum.push_back( s.getWeight(step) );
        std::cout << step << " MeV: " << _spectrum.back() << std::endl;
      }
    }

  private:

    std::vector<double> _spectrum;

  };

} // end of namespace mu2e

#endif /* Mu2eUtilities_BinnedSpectrum_hh */
