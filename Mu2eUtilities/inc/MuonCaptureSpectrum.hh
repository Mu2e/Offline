#ifndef Mu2eUtilities_MuonCaptureSpectrum_hh
#define Mu2eUtilities_MuonCaptureSpectrum_hh

// Mu2e includes
#include "Mu2eUtilities/inc/Table.hh"

// CLHEP includes
#include "CLHEP/Vector/LorentzVector.h"

// C++ includes
#include <utility>

namespace mu2e {

  class MuonCaptureSpectrum {

  public:
    
    enum enum_type    { Flat  , RMC      };
    enum enum_type_2D { Flat2D, KrollWadaJoseph };


    MuonCaptureSpectrum() : _spectrum( RMC ), _spectrum2D( KrollWadaJoseph ) {}
    ~MuonCaptureSpectrum(){}
    MuonCaptureSpectrum(const bool kMaxUserSet, const double kMaxUser, const double kMaxMax) : _kMaxUserSet( kMaxUserSet), _kMaxUser( kMaxUser), _kMaxMax ( kMaxMax) {}

    double getWeight   ( const double E ) const;
    double get2DWeight ( const double x, const double y, const double E ) const; 
    double get2DMax    ( const double E ) const;


 
    void   setSpectrum   ( enum_type    spectrum   ) { _spectrum   = spectrum;   }
    void   setSpectrum2D ( enum_type_2D spectrum2D ) { _spectrum2D = spectrum2D; }


    static double getFlat( const double e, const double x = 0., const double y = 0.);
    static double getRMCSpectrum( const double e , const bool kMaxUserSet, const double kMaxUser, const double kMaxMax);
    static std::pair<CLHEP::HepLorentzVector,CLHEP::HepLorentzVector> getElecPosiVectors( const double energy,
                                                                                          const double x,
                                                                                          const double y );

    double getKrollWadaJosephSpectrum( const double e, const double x, const double y ) const;


  private:

    enum_type    _spectrum{};
    enum_type_2D _spectrum2D{};
    bool _kMaxUserSet{false};
    double _kMaxUser{0.};
    double _kMaxMax{0};

  };

} // end of namespace mu2e

#endif /* Mu2eUtilities_MuonCaptureSpectrum_hh */

