#ifndef Mu2eUtilities_SimpleSpectrum_hh
#define Mu2eUtilities_SimpleSpectrum_hh
//
// Original Author: Kyle Knoepfel
//

// Read Simple DIO spectrum from a table and merge it
// with the spectrum coming from the endpoint region formula
//
// Three spectrum options are available:
//
//   (1) "flat" : where electrons are produced uniformly in energy, up
//                to the endpoint energy
//   (2) "pol5" : the order-5 polynomial approximation given by
//                Czarnecki, et al (PRD 84, 013006)
//   (3) "pol58": the order-5 through 8 polynomial approximation given
//                by Czarnecki, et al (PRD 84, 013006)

// C++ includes
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "GeneralUtilities/inc/EnumToStringSparse.hh"
// Mu2e includes
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/PhysicsParams.hh"

namespace mu2e {

  class SimpleSpectrum {

  public:
    
    class SpectrumType {
    public:
      enum enum_type { unknown, Flat, FlatTrunc, Pol5, Pol58 };
      static std::string const& typeName() {
        static std::string type("SimpleSpectrumType"); return type;
      }
      static std::map<enum_type,std::string> const& names() {
        static std::map<enum_type,std::string> nam;
        
        if ( nam.empty() ) {
          nam[unknown]     = "unknown";
          nam[Flat]        = "flat";
          nam[FlatTrunc]   = "flatTrunc";
          nam[Pol5]        = "pol5";
          nam[Pol58]       = "pol58";
        }

        return nam;
      }
    };

    typedef EnumToStringSparse<SpectrumType> Spectrum;

    SimpleSpectrum( Spectrum::enum_type i, const PhysicsParams& phy = *GlobalConstantsHandle<PhysicsParams>() ) ;
    ~SimpleSpectrum();

    double getWeight(double e) const;

    static double getFlat     (const double e, const PhysicsParams& phy = *GlobalConstantsHandle<PhysicsParams>() );
    static double getFlatTrunc(const double e, const PhysicsParams& phy = *GlobalConstantsHandle<PhysicsParams>() );
    static double getPol5     (const double e, const PhysicsParams& phy = *GlobalConstantsHandle<PhysicsParams>() );
    static double getPol58    (const double e, const PhysicsParams& phy = *GlobalConstantsHandle<PhysicsParams>() );

  private:

    Spectrum::enum_type _approx;
    PhysicsParams _phy;

  };

} // end of namespace mu2e

#endif /* Mu2eUtilities_SimpleSpectrum_hh */

