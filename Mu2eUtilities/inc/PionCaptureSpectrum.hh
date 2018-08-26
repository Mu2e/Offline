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
// #include "Mu2eUtilities/inc/Table.hh"

// CLHEP includes
// #include "CLHEP/Vector/LorentzVector.h"

// C++ includes
#include <utility>

namespace CLHEP {
  class RandFlat;
  class HepLorentzVector;
}

namespace mu2e {

  class RandomUnitSphere;

  class PionCaptureSpectrum {
  public:

    enum enum_type    { Flat  , Bistirlich      };
    enum enum_type_2D { Flat2D, KrollWadaJoseph };

    PionCaptureSpectrum() : _spectrum( Bistirlich ), _spectrum2D( KrollWadaJoseph ) {}
    ~PionCaptureSpectrum(){}

    double getWeight   (double E) const;
    double get2DWeight (double x, double y, double E) const;
    double get2DMax    (double E) const;

    void   setSpectrum   ( enum_type    spectrum   ) { _spectrum   = spectrum;   }
    void   setSpectrum2D ( enum_type_2D spectrum2D ) { _spectrum2D = spectrum2D; }

    double getFlat              (double e, double x = 0., double y = 0.) const ;
    double getBistirlichSpectrum(double e) const ;

    void   getElecPosiVectors(RandomUnitSphere* randomUnitSphere,
			      CLHEP::RandFlat*  randFlat,
			      double energy,
			      double x,
			      double y,
			      CLHEP::HepLorentzVector& mome,
			      CLHEP::HepLorentzVector& momp) const;

    double getKrollWadaJosephSpectrum( const double e, const double x, const double y ) const;

  private:

    enum_type    _spectrum;
    enum_type_2D _spectrum2D;

  };

} // end of namespace mu2e

#endif /* Mu2eUtilities_PionCaptureSpectrum_hh */
