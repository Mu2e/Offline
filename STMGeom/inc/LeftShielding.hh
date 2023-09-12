#ifndef STMGeom_LeftShielding_hh
#define STMGeom_LeftShielding_hh

// Germanium Detector Object
//
// Author: Haichuan Cao
// Sept 2023

#include <string>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class LeftShielding {
  public:
    LeftShielding(bool build, double L_Length,
    double Pb_depth, double Cu_depth, double BP_depth, double Left_Xmin,
    CLHEP::Hep3Vector const & originInMu2e = CLHEP::Hep3Vector(), CLHEP::HepRotation const & rotation = CLHEP::HepRotation()
   ):
      _build(build), _L_Length(L_Length),
      _Lleaddepth(Pb_depth), _Lcopperdepth(Cu_depth), _LBPdepth(BP_depth), _Left_Xmin(Left_Xmin),
      _originInMu2e(originInMu2e), _rotation(rotation)
    {
    }

   bool    build()            const {return _build;}
   double  L_Length()         const {return _L_Length;}
   double  Lleaddepth()       const {return _Lleaddepth;}
   double  Lcopperdepth()     const {return _Lcopperdepth;}
   double  LBPdepth()         const {return _LBPdepth;}
   double  Left_Xmin()        const {return _Left_Xmin;}


   double  Left_Thickness()   const {return _Lleaddepth*2 + _LBPdepth*2 + _Lcopperdepth;}

   LeftShielding() {}

  private:

    bool               _build;
    double             _L_Length;
    double             _Lleaddepth;
    double             _Lcopperdepth;
    double             _LBPdepth;
    double             _Left_Xmin;

    CLHEP::Hep3Vector  _originInMu2e;
    CLHEP::HepRotation _rotation;
  };

}

#endif/*STMGeom_LeftShielding_hh*/
