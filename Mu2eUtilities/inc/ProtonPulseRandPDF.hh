#ifndef ProtonPulseRandPDF_hh
#define ProtonPulseRandPDF_hh

//
// Constructor of a PDF to extract random times to describe the proton pulse
//
//
// Original author: Gianni Onorato
//                  Kyle Knoepfel (significant updates)
//

// Mu2e includes
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "Offline/GeneralUtilities/inc/EnumToStringSparse.hh"
#include "Offline/Mu2eUtilities/inc/Table.hh"

// CLHEP includes
#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Random/RandomEngine.h"

// Framework includes
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
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

    // Required for EnumToStringSparse
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

    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::Atom<std::string> pulseType {
        Name("pulseType"),
          Comment("Allowed values are: default, total, oot, allflat.")
          };

      fhicl::Atom<double> limitingHalfWidth {
        Name("limitingHalfWidth"),
          Comment("limit on the shape tail during generation")
          };

      fhicl::Atom<std::string> potPulse {
        Name("potPulse"),
          Comment("Text file with proton pulse shape"),
          };

      fhicl::Atom<std::string> acDipole {
        Name("acDipole"),
          Comment("Text file with AC dipole transmission function"),
          };

      fhicl::OptionalAtom<double> tmin {
        Name("tmin"),
          Comment("Override proton pulse start time.  The default is determined from the ConditionsService based on the pulseType.")
          };

      fhicl::OptionalAtom<double> tmax {
        Name("tmax"),
          Comment("Override proton pulse end time.  The default is determined from the ConditionsService based on the pulseType.")
          };

      fhicl::Atom<double> tres {
        Name("tres"),
          Comment("The time bin size, in ns."),
          1.
          };
    };

    ProtonPulseRandPDF(art::RandomNumberGenerator::base_engine_t& engine, const Config& conf);

    ~ProtonPulseRandPDF(){}

    std::string pulseType()                  const { return pulseEnum_.name(); }
    double fire();
    const std::vector<double>& getSpectrum() const { return spectrum_; }
    const std::vector<double>& getTimes()    const { return times_   ; }

  private:

    const ProtonPulseEnum pulseEnum_;

    double limitingHalfWidth_;
    double DRPeriod_;

    double tmin_;
    double tmax_;
    double tres_;
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
