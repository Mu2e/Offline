#ifndef ProtonPulseRandPDF_hh
#define ProtonPulseRandPDF_hh

//  
// Constructor of a PDF to extract random times to describe the proton pulse
//
// $Id: ProtonPulseRandPDF.hh,v 1.11 2014/04/25 17:26:42 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2014/04/25 17:26:42 $
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
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

// C++ includes
#include <map>
#include <vector>

namespace mu2e {

  // ===================== Enum detail base class =====================
  class ProtonPulseEnumDetail {
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
  
  typedef EnumToStringSparse<ProtonPulseEnumDetail> ProtonPulseEnum;

  // ==================================================================    
  // ProtonPulseRandPDF class declaration
  // ==================================================================    
  class ProtonPulseRandPDF : public ProtonPulseEnum {
  public:
    ProtonPulseRandPDF(art::RandomNumberGenerator::base_engine_t& engine,
                       const fhicl::ParameterSet pset );
    ~ProtonPulseRandPDF(){}

    std::string pulseType()                  const { return pulseEnum_.name(); }
    double fire();
    const std::vector<double>& getSpectrum() const { return spectrum_; }
    const std::vector<double>& getTimes()    const { return times_   ; }

  private:

    const AcceleratorParams* accPar_;

    const ProtonPulseEnum pulseEnum_;
    
    const double tmin_;
    const double tmax_;
    const double tres_;
    const std::vector<double> times_;
    
    const double extFactor_;    
    const Table<2> pulseShape_;
    const Table<2> acdipole_;

    const std::vector<double> spectrum_;

    CLHEP::RandGeneral randSpectrum_;

    // Modifiers
    double setTmin();
    double setTmax();
    std::vector<double> setTimes();
    std::vector<double> setSpectrum();

    double determineIntrinsicExt( const std::string& filename );
    TableVec<2> setPotPulseShape( const std::string& filename );

  };

}

#endif /* ProtonPulseRandPDF_hh */
