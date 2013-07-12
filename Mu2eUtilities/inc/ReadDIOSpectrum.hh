#ifndef Mu2eUtilities_ReadDIOSpectrum_hh
#define Mu2eUtilities_ReadDIOSpectrum_hh
//
// Generate a momentum for the DIO electrons, using custom
// spectrum 
//
// $Id: ReadDIOSpectrum.hh,v 1.4 2013/07/12 17:17:38 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2013/07/12 17:17:38 $
//
// Original Author: Gianni Onorato
//                  Kyle Knoepfel (significant updates)

// Mu2e includes
#include "Mu2eUtilities/inc/DIOBase.hh"

// C++ includes
#include <memory>
#include <vector>

namespace mu2e {

  class DIOBase;

  class ReadDIOSpectrum {

  public:

    ReadDIOSpectrum( double emin, double emax, double spectRes, std::string spectrum ); 
    ~ReadDIOSpectrum();

    // To make accessible to CLHEP::RandGeneral, return pointer to 
    // first spectrum entry
    const double* getPDF()   const { return &(*_spectrum.begin() ); }
    unsigned      getNbins() const { return _spectrum.size();       }

  private:

    double _emin; 
    double _emax; 
    double _res;
    std::string _spectrumName;

    std::unique_ptr<DIOBase> _dioShape;
    std::vector<double> _spectrum;

    inline void makeSpectrum() {
      for ( double step = _emin ; step <= _emax ; step += _res ) 
        _spectrum.push_back( _dioShape->getWeight(step) );
    }

  };

} // end of namespace mu2e

#endif /* Mu2eUtilities_ReadDIOSpectrum_hh */

