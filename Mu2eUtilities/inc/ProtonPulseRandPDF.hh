#ifndef ProtonPulseRandPDF_hh
#define ProtonPulseRandPDF_hh

//  
// Constructor of a PDF to extract random times to describe the proton pulse
//
// $Id: ProtonPulseRandPDF.hh,v 1.5 2014/02/07 14:48:44 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2014/02/07 14:48:44 $
//
// Original author: Gianni Onorato
//                  Kyle Knoepfel (significant updates)
//

// Mu2e includes
#include "Mu2eUtilities/inc/Table.hh"

// CLHEP includes
#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Random/RandomEngine.h"

// Framework includes
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

// C++ includes
#include <vector>

namespace mu2e {

  class ProtonPulseRandPDF {

  public:

    ProtonPulseRandPDF(art::RandomNumberGenerator::base_engine_t& engine);
    ~ProtonPulseRandPDF(){}

    double fire();
    const std::vector<double>& getSpectrum() const { return _spectrum; }

  private:

    const Table<2> _pulseShape;
    const Table<2> _acdipole;
    std::vector<double> _spectrum;
    double _timeMin;
    double _timeMax;
    CLHEP::RandGeneral _randSpectrum;

    //PDF description
    std::vector<double> setSpectrum();

  };

}

#endif /* ProtonPulseRandPDF_hh */
