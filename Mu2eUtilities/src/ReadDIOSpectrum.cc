//
// Generate an energy value for the DIO electrons, using a generic 
// custom distribution as momentum spectrum
// The construction
// of the spectrum is made by specialized classes
//
// $Id: ReadDIOSpectrum.cc,v 1.5 2013/07/12 17:17:38 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2013/07/12 17:17:38 $
//
// Original author: Gianni Onorato
//                  Kyle Knoepfel (significant updates)

// C++ includes
#include <iostream>

// Framework includes
#include "cetlib/exception.h"

// Mu2e includes
#include "Mu2eUtilities/inc/ReadDIOSpectrum.hh"
#include "Mu2eUtilities/inc/CzarneckiSpectrum.hh"
#include "Mu2eUtilities/inc/ShankerWatanabeSpectrum.hh"
#include "Mu2eUtilities/inc/SimpleSpectrum.hh"

using namespace std;

namespace mu2e {

  ReadDIOSpectrum::ReadDIOSpectrum(double emin, double emax, double spectRes, string spectrum ) :
    _emin ( emin ),
    _emax ( emax ),
    _res ( spectRes ),
    _spectrumName( spectrum ) 
  {

    if      ( _spectrumName == "ShankerWatanabe" ) _dioShape.reset( new ShankerWatanabeSpectrum() );
    else if ( _spectrumName == "Czarnecki"       ) _dioShape.reset( new CzarneckiSpectrum()       );
    else if ( _spectrumName == "flat"            ) _dioShape.reset( new SimpleSpectrum( SimpleSpectrum::Spectrum::Flat  )  );
    else if ( _spectrumName == "pol5"            ) _dioShape.reset( new SimpleSpectrum( SimpleSpectrum::Spectrum::Pol5  )  );
    else if ( _spectrumName == "pol58"           ) _dioShape.reset( new SimpleSpectrum( SimpleSpectrum::Spectrum::Pol58 )  );
    else {
      throw cet::exception("MODEL")
        << "Wrong or not allowed DIO energy spectrum";
    }

    makeSpectrum();

  }

  ReadDIOSpectrum::~ReadDIOSpectrum() {
  }

}

