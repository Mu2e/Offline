#ifndef ProtonPulseRandPDF_hh
#define ProtonPulseRandPDF_hh

//  
// Constructor of a PDF to extract random times to describe the proton pulse
//
// $Id: ProtonPulseRandPDF.hh,v 1.9 2014/04/01 15:03:16 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2014/04/01 15:03:16 $
//
// Original author: Gianni Onorato
//                  Kyle Knoepfel (significant updates)
//

// Mu2e includes
#include "ConditionsService/inc/AcceleratorParams.hh"
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
      enum enum_type { unknown, DEFAULT, TOTAL, OOT, OOTFLAT, ALLFLAT };
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
          nam[OOTFLAT] = "ootflat";
          nam[ALLFLAT] = "allflat";
        }

        return nam;
      }
    };

    typedef EnumToStringSparse<PotSpectrumType> PotSpectrum;
    
    ProtonPulseRandPDF(art::RandomNumberGenerator::base_engine_t& engine,
                       const std::string pulseString = "default" );
    ~ProtonPulseRandPDF(){}

    double fire();
    const std::vector<double>& getSpectrum() const { return spectrum_; }
    const std::vector<double>& getTimes()    const { return times_; }

  private:

    const AcceleratorParams* accPar_;

    std::vector<double> potSpectrum_;
    std::vector<double> dipoleSpectrum_;
    std::vector<double> ootSpectrum_;

    std::vector<double> times_;

    double timeMin_;
    double timeMax_;
    double fireOffset_;

    const Table<2> pulseShape_;
    const Table<2> acdipole_;
    const Table<2> ootPulse_;

    const PotSpectrum pulseEnum_;
    const std::size_t nPoints_;

    const double extFactor_;

    std::vector<double> spectrum_;

    CLHEP::RandGeneral randSpectrum_;

    std::size_t calculateNpoints();

    //PDF description
    std::vector<double> setSpectrum() const;
    
    std::vector<double> getShape( const Table<2>& table, const double timeOffset = 0. ) const;
    void renormalizeShape( std::vector<double>& shape, const double norm ) const;

  };

}

#endif /* ProtonPulseRandPDF_hh */
