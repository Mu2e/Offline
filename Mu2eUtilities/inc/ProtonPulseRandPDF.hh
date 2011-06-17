#ifndef ProtonPulseRandPDF_hh
#define ProtonPulseRandPDF_hh

//  
// Constructor of a PDF to extract random times to describe the proton pulse
//
// $Id: ProtonPulseRandPDF.hh,v 1.3 2011/06/17 21:08:27 onoratog Exp $
// $Author: onoratog $
// $Date: 2011/06/17 21:08:27 $
//
// Original author Gianni Onorato
//
//

#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Random/RandomEngine.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include <vector>

namespace mu2e {

  class ProtonPulseRandPDF {

  public:

    ProtonPulseRandPDF(art::RandomNumberGenerator::base_engine_t& engine);

    ~ProtonPulseRandPDF();

    double fire();

    int calculateNBins();

  private:

    int _nBinsSpectrum;

    CLHEP:: RandGeneral _randSpectrum;

    //PDF description
    std::vector<double> ProtonPulseSpectrum();


  };

}

#endif /* ProtonPulseRandPDF_hh */
