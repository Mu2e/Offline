#ifndef Mu2eUtilities_ReadDIOSpectrum_hh
#define Mu2eUtilities_ReadDIOSpectrum_hh
//
// Generate a momentum for the DIO electrons, using custom
// spectrum (so far Shanker-Wanatabe or Czarnecki
//
// $Id: ReadDIOSpectrum.hh,v 1.1 2012/02/06 23:56:32 onoratog Exp $
// $Author: onoratog $
// $Date: 2012/02/06 23:56:32 $
//
// Original Author: Gianni Onorato
//

// C++ includes
#include <vector>

// Framework includes
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

// Mu2e includes
#include "Mu2eUtilities/inc/DIOBase.hh"

//CLHEP includes
#include "CLHEP/Random/RandGeneral.h"

namespace mu2e {

  class ReadDIOSpectrum: public DIOBase {

  public:

    ReadDIOSpectrum(int atomicZ, double emin, double emax, double spectRes,
		    std::string spectrum, art::RandomNumberGenerator::base_engine_t & engine);

    ~ReadDIOSpectrum();

    double fire();

  private:

    int _Znum;

    double _emin, _emax, _res;

    int calculateNBins();

    int _nBinsSpectrum;

    std::string _spectrum;

    CLHEP::RandGeneral _randEnergy;

    std::vector<double> ReadSpectrum();

  };

} // end of namespace mu2e

#endif /* Mu2eUtilities_ReadDIOSpectrum_hh */

