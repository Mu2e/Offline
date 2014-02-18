#ifndef ProtonPulseRandPDF_hh
#define ProtonPulseRandPDF_hh

//  
// Constructor of a PDF to extract random times to describe the proton pulse
//
// $Id: ProtonPulseRandPDF.hh,v 1.6 2014/02/18 20:21:25 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2014/02/18 20:21:25 $
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

    enum enum_type { DEFAULT, TOTAL, OOT };

    ProtonPulseRandPDF(art::RandomNumberGenerator::base_engine_t& engine,
                       const std::string pulseString = "default" );
    ~ProtonPulseRandPDF(){}

    double fire();
    const std::vector<double>& getSpectrum() const { return _spectrum; }

  private:

    std::vector<double> _potSpectrum;
    std::vector<double> _dipoleSpectrum;
    std::vector<double> _ootSpectrum;

    std::vector<double> _times;

    double _timeMin;
    double _timeMax;
    double _fireOffset;

    const Table<2> _pulseShape;
    const Table<2> _acdipole;
    const Table<2> _ootPulse;

    const enum_type _pulseEnum;
    const std::size_t _nPoints;

    const double _extFactor;

    std::vector<double> _spectrum;

    CLHEP::RandGeneral _randSpectrum;

    enum_type getPulseEnum( const std::string& pulseString ) const;

    std::size_t calculateNpoints();

    //PDF description
    std::vector<double> setSpectrum() const;
    
    std::vector<double> getShape( const Table<2>& table, const double timeOffset = 0. ) const;
    void renormalizeShape( std::vector<double>& shape, const double norm ) const;

  };

}

#endif /* ProtonPulseRandPDF_hh */
