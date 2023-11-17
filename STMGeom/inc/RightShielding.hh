#ifndef STMGeom_RightShielding_hh
#define STMGeom_RightShielding_hh

// Germanium Detector Object
//
// Author: Haichuan Cao
// Sept 2023

#include <string>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class RightShielding {
  public:
    RightShielding(bool build, double R_Length,
    double Pb_depth, double Cu_depth, double BP_depth, double Right_Xmax,
    CLHEP::Hep3Vector const & originInMu2e = CLHEP::Hep3Vector(), CLHEP::HepRotation const & rotation = CLHEP::HepRotation()
   ):
      _build(build), _R_Length(R_Length),
      _Rleaddepth(Pb_depth), _Rcopperdepth(Cu_depth), _RBPdepth(BP_depth), _Right_Xmax(Right_Xmax),
      _originInMu2e(originInMu2e), _rotation(rotation)
    {
    }

   bool    build()            const {return _build;}
   double  R_Length()         const {return _R_Length;}
   double  Rleaddepth()       const {return _Rleaddepth;}
   double  Rcopperdepth()     const {return _Rcopperdepth;}
   double  RBPdepth()         const {return _RBPdepth;}
   double  Right_Xmax()       const {return _Right_Xmax;}


   double  Right_Thickness()   const {return (_Rleaddepth*2 + _RBPdepth*2 + _Rcopperdepth)*sqrt(2);}

   RightShielding() {}

  private:

    bool               _build;
    double             _R_Length;
    double             _Rleaddepth;
    double             _Rcopperdepth;
    double             _RBPdepth;
    double             _Right_Xmax;

    CLHEP::Hep3Vector  _originInMu2e;
    CLHEP::HepRotation _rotation;
  };

}

#endif/*STMGeom_RightShielding_hh*/
