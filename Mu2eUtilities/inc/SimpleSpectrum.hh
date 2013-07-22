#ifndef Mu2eUtilities_SimpleSpectrum_hh
#define Mu2eUtilities_SimpleSpectrum_hh
//
// Read Simple DIO spectrum from a table and merge it
// with the spectrum coming from the endopoint region formula

// $Id: SimpleSpectrum.hh,v 1.2 2013/07/22 18:57:42 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2013/07/22 18:57:42 $
//
// Original Author: Kyle Knoepfel
//

// C++ includes
#include <map>
#include <utility>
#include <vector>

// Mu2e includes
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

    SimpleSpectrum( Spectrum::enum_type i) ;
    ~SimpleSpectrum();

    double getWeight(double e);

    static double getFlat (double e);
    static double getPol5 (double e);
    static double getPol58(double e);

  private:

    Spectrum::enum_type _approx;

  };

} // end of namespace mu2e

#endif /* Mu2eUtilities_SimpleSpectrum_hh */

