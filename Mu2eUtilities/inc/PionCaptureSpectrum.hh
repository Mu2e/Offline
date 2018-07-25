#ifndef Mu2eUtilities_PionCaptureSpectrum_hh
#define Mu2eUtilities_PionCaptureSpectrum_hh
//
// Read PionCapture DIO spectrum from a table and merge it
// with the spectrum coming from the endopoint region formula

// $Id: PionCaptureSpectrum.hh,v 1.1 2014/02/05 21:32:26 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2014/02/05 21:32:26 $
//
// Original Author: Kyle Knoepfel
//

// Mu2e includes
#include "Mu2eUtilities/inc/Table.hh"

// CLHEP includes
#include "CLHEP/Vector/LorentzVector.h"

// C++ includes
#include <utility>

namespace CLHEP {
  class RandFlat;
}

namespace mu2e {

  class RandomUnitSphere;

  class PionCaptureSpectrum {
  public:

    enum enum_type    { Flat  , Bistirlich      };
    enum enum_type_2D { Flat2D, KrollWadaJoseph };

    PionCaptureSpectrum() : _spectrum( Bistirlich ), _spectrum2D( KrollWadaJoseph ) {}
    ~PionCaptureSpectrum(){}

    double getWeight   ( const double E ) const;
    double get2DWeight ( const double x, const double y, const double E ) const;
    double get2DMax    ( const double E ) const;

    void   setSpectrum   ( enum_type    spectrum   ) { _spectrum   = spectrum;   }
    void   setSpectrum2D ( enum_type_2D spectrum2D ) { _spectrum2D = spectrum2D; }

    static double getFlat( const double e, const double x = 0., const double y = 0.);
    static double getBistirlichSpectrum( const double e );
    static std::pair<CLHEP::HepLorentzVector,CLHEP::HepLorentzVector>
    getElecPosiVectors(RandomUnitSphere& randomUnitSphere,
                       CLHEP::RandFlat& randFlat,
                       const double energy,
                       const double x,
                       const double y);

    double getKrollWadaJosephSpectrum( const double e, const double x, const double y ) const;

  private:

    enum_type    _spectrum;
    enum_type_2D _spectrum2D;

  };

} // end of namespace mu2e

#endif /* Mu2eUtilities_PionCaptureSpectrum_hh */
