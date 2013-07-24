#ifndef Mu2eUtilities_SimpleSpectrum_hh
#define Mu2eUtilities_SimpleSpectrum_hh
//
// Read Simple DIO spectrum from a table and merge it
// with the spectrum coming from the endopoint region formula

// $Id: SimpleSpectrum.hh,v 1.3 2013/07/24 18:48:24 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2013/07/24 18:48:24 $
//
// Original Author: Kyle Knoepfel
//

// C++ includes
#include <map>
#include <utility>
#include <vector>

// Mu2e includes
#include "ConditionsService/inc/GlobalConstantsHandle.hh"
#include "ConditionsService/inc/PhysicsParams.hh"
#include "GeneralUtilities/inc/EnumToStringSparse.hh"

namespace mu2e {

  class SimpleSpectrum {

  public:
    
    class SpectrumType {
    public:
      enum enum_type { unknown, Flat, Pol5, Pol58 };
      static std::string const& typeName() {
        static std::string type("SimpleSpectrumType"); return type;
      }
      static std::map<enum_type,std::string> const& names() {
        static std::map<enum_type,std::string> nam;
        
        if ( nam.empty() ) {
          nam[unknown] = "unknown";
          nam[Flat]    = "flat";
          nam[Pol5]    = "pol5";
          nam[Pol58]   = "pol58";
        }

        return nam;
      }
    };

    typedef EnumToStringSparse<SpectrumType> Spectrum;

    SimpleSpectrum( Spectrum::enum_type i, const PhysicsParams& phy = *GlobalConstantsHandle<PhysicsParams>() ) ;
    ~SimpleSpectrum();

    double getWeight(double e);

    static double getFlat (const double e, const PhysicsParams& phy = *GlobalConstantsHandle<PhysicsParams>() );
    static double getPol5 (const double e, const PhysicsParams& phy = *GlobalConstantsHandle<PhysicsParams>() );
    static double getPol58(const double e, const PhysicsParams& phy = *GlobalConstantsHandle<PhysicsParams>() );

  private:

    Spectrum::enum_type _approx;
    PhysicsParams _phy;

  };

} // end of namespace mu2e

#endif /* Mu2eUtilities_SimpleSpectrum_hh */

