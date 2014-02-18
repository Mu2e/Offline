#ifndef ProtonPulseRandPDF_hh
#define ProtonPulseRandPDF_hh

//  
// Constructor of a PDF to extract random times to describe the proton pulse
//
// $Id: ProtonPulseRandPDF.hh,v 1.7 2014/02/18 22:08:44 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2014/02/18 22:08:44 $
//
// Original author: Gianni Onorato
//                  Kyle Knoepfel (significant updates)
//

// Mu2e includes
#include "GeneralUtilities/inc/EnumToStringSparse.hh"
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

    class PotSpectrumType {
    public:
      enum enum_type { unknown, DEFAULT, TOTAL, OOT };
      static std::string const& typeName() {
        static std::string type("PotSpectrumType"); return type;
      }
      static std::map<enum_type,std::string> const& names() {
        static std::map<enum_type,std::string> nam;
        
        if ( nam.empty() ) {
          nam[unknown] = "unknown";
          nam[DEFAULT] = "default";
          nam[TOTAL]   = "total";
          nam[OOT]     = "oot";
        }

        return nam;
      }
    };

    typedef EnumToStringSparse<PotSpectrumType> PotSpectrum;
    
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

    const PotSpectrum _pulseEnum;
    const std::size_t _nPoints;

    const double _extFactor;

    std::vector<double> _spectrum;

    CLHEP::RandGeneral _randSpectrum;

    std::size_t calculateNpoints();

    //PDF description
    std::vector<double> setSpectrum() const;
    
    std::vector<double> getShape( const Table<2>& table, const double timeOffset = 0. ) const;
    void renormalizeShape( std::vector<double>& shape, const double norm ) const;

  };

}

#endif /* ProtonPulseRandPDF_hh */
