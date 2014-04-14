#ifndef ProtonPulseRandPDF_hh
#define ProtonPulseRandPDF_hh

//  
// Constructor of a PDF to extract random times to describe the proton pulse
//
// $Id: ProtonPulseRandPDF.hh,v 1.10 2014/04/14 18:12:55 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2014/04/14 18:12:55 $
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
#include <map>
#include <vector>

namespace mu2e {

  class ProtonPulseRandPDF {

  public:

    class PotSpectrumType {
    public:
      enum enum_type { unknown, DEFAULT, TOTAL, OOT, ALLFLAT };
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
          nam[ALLFLAT] = "allflat";
        }

        return nam;
      }
    };

    typedef EnumToStringSparse<PotSpectrumType> PotSpectrumEnum;
    
    ProtonPulseRandPDF(art::RandomNumberGenerator::base_engine_t& engine,
                       const std::string pulseString = "default" );
    ~ProtonPulseRandPDF(){}

    double fire();
    const std::vector<double>& getSpectrum() const { return spectrum_; }
    const std::vector<double>& getTimes()    const { return times_   ; }

  private:

    const AcceleratorParams* accPar_;

    std::map<double,double> potSpectrum_;
    std::map<double,double> dipoleSpectrum_;

    double timeMin_;
    double timeMax_;
    double fireOffset_;

    const Table<2> pulseShape_;
    const Table<2> acdipole_;
    
    const double extFactor_;
    
    const PotSpectrumEnum pulseEnum_;
    
    std::vector<double> times_;
    std::vector<double> spectrum_;

    CLHEP::RandGeneral randSpectrum_;

    std::vector<double> setTimes();

    //PDF description
    std::vector<double> setSpectrum() const;

    std::vector<double> preparePotSpectrum() const;

    std::map<double,double> getShape   ( const Table<2>& table, const double timeOffset = 0. ) const;

    static void   renormalizeShape     ( std::map<double,double>& shape      , const double norm );
    static double determineIntrinsicExt( const std::map<double,double>& shape, const double hw   );
    static void   replaceOotShape      ( std::map<double,double>& shape      , const double hw, const double norm );

  };

}

#endif /* ProtonPulseRandPDF_hh */
