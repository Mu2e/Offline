#ifndef Mu2eUtilities_MuonCaptureSpectrum_hh
#define Mu2eUtilities_MuonCaptureSpectrum_hh

// Mu2e includes
// #include "Mu2eUtilities/inc/Table.hh"

// CLHEP includes
#include "CLHEP/Vector/LorentzVector.h"

// C++ includes
#include <utility>


namespace CLHEP {
  class RandFlat;
}
namespace mu2e {

  class RandomUnitSphere;

  class MuonCaptureSpectrum {

  public:
    
    enum enum_type    { Flat  , RMC             };
    enum enum_type_2D { Flat2D, KrollWadaJoseph };

   MuonCaptureSpectrum(){}

    // random number generators ar owned by the callers, no memory cleanup needed
    MuonCaptureSpectrum(CLHEP::RandFlat* randFlat, RandomUnitSphere* randomUnitSphere);

    MuonCaptureSpectrum(bool kMaxUserSet, double kMaxUser, double kMaxMax,
			CLHEP::RandFlat* randFlat = 0, RandomUnitSphere* randomUnitSphere = 0);

    ~MuonCaptureSpectrum(){}

    double getWeight   (double E) const;
    double get2DWeight (double x, double y, double E) const; 
    double get2DMax    (double E) const;
 
    void   setSpectrum   (enum_type    spectrum  ) { _spectrum   = spectrum;   }
    void   setSpectrum2D (enum_type_2D spectrum2D) { _spectrum2D = spectrum2D; }


    double getFlat       (double e, double x = 0., double y = 0.) const ;
    double getRMCSpectrum(double e, bool kMaxUserSet, double kMaxUser, double kMaxMax) const;

    void   getElecPosiVectors(double energy, CLHEP::HepLorentzVector& mome, CLHEP::HepLorentzVector& momp) const;

    double getKrollWadaJosephSpectrum(double e, double x, double y) const;

    void   fire(double energy, double& x, double& y) const; 
   

  private:

    enum_type          _spectrum;
    enum_type_2D       _spectrum2D;
    bool               _kMaxUserSet;
    double             _kMaxUser;
    double             _kMaxMax;

    CLHEP::RandFlat*   _rnFlat;
    RandomUnitSphere*  _rnUnitSphere;

    double             _me;		// electron mass
    double             _mmu;		// muon mass
    double             _MN;		// mass of the initial state nucleus
  };

} // end of namespace mu2e

#endif /* Mu2eUtilities_MuonCaptureSpectrum_hh */

