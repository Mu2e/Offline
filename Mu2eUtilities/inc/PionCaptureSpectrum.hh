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

    // random number generators ar owned by the callers, no memory cleanup needed
    PionCaptureSpectrum(CLHEP::RandFlat* randFlat, RandomUnitSphere* randomUnitSphere);

    PionCaptureSpectrum(bool kMaxUserSet, double kMaxUser, double kMaxMax,
			CLHEP::RandFlat* randFlat, RandomUnitSphere* randomUnitSphere);

    ~PionCaptureSpectrum(){}

    double getWeight   (double E) const;
    double get2DWeight (double x, double y, double E) const;
    double get2DMax    (double E) const;

    void   setSpectrum   ( enum_type    spectrum   ) { _spectrum   = spectrum;   }
    void   setSpectrum2D ( enum_type_2D spectrum2D ) { _spectrum2D = spectrum2D; }

    double getFlat              (double e, double x = 0., double y = 0.) const ;
    double getBistirlichSpectrum(double e) const ;

    void   fire(double energy, double& x, double& y) const; 

    void   getElecPosiVectors(double energy, CLHEP::HepLorentzVector& mome, CLHEP::HepLorentzVector& momp) const;

    double getKrollWadaJosephSpectrum( const double e, const double x, const double y ) const;

  private:

    enum_type          _spectrum;
    enum_type_2D       _spectrum2D;
    bool               _kMaxUserSet;
    double             _kMaxUser;
    double             _kMaxMax;

    CLHEP::RandFlat*   _rnFlat;
    RandomUnitSphere*  _rnUnitSphere;

    double             _me;		// electron mass
    double             _mpi;		// pi- mass
    double             _MN;		// mass of the initial state nucleus (Al)
  };

} // end of namespace mu2e

#endif /* Mu2eUtilities_PionCaptureSpectrum_hh */
