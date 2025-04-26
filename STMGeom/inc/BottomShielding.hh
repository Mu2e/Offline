#ifndef STMGeom_BottomShielding_hh
#define STMGeom_BottomShielding_hh

// Germanium Detector Object
//
// Author: Haichuan Cao
// Sept 2023

#include <string>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class BottomShielding {
  public:
    BottomShielding(bool build, double floor_Zlength, double Front_LB, double Front_LB_inner,
    double Pb_depth, double Cu_depth, double BP_depth,
    CLHEP::Hep3Vector const & originInMu2e = CLHEP::Hep3Vector(), CLHEP::HepRotation const & rotation = CLHEP::HepRotation()
   ):
      _build(build), _floor_Zlength(floor_Zlength), _Front_LB(Front_LB), _Front_LB_inner(Front_LB_inner),
      _Bleaddepth(Pb_depth), _Bcopperdepth(Cu_depth), _BBPdepth(BP_depth),
      _originInMu2e(originInMu2e), _rotation(rotation)
    {
    }

   bool    build()            const {return _build;}
   double  floor_Zlength()    const {return _floor_Zlength;}
   double  Front_LB()         const {return _Front_LB;}
   double  Front_LB_inner()   const {return _Front_LB_inner;}

   double  Bleaddepth()       const {return _Bleaddepth;}
   double  Bcopperdepth()     const {return _Bcopperdepth;}
   double  BBPdepth()         const {return _BBPdepth;}

   double  Bottom_Thickness()  const {return _Bleaddepth*2 + _BBPdepth*2 + _Bcopperdepth;}

   BottomShielding() {}

  private:

    bool               _build;
    double             _floor_Zlength;
    double             _Front_LB;
    double             _Front_LB_inner;
    double             _Bleaddepth;
    double             _Bcopperdepth;
    double             _BBPdepth;

    CLHEP::Hep3Vector  _originInMu2e;
    CLHEP::HepRotation _rotation;
  };

}

#endif/*STMGeom_BottomShielding_hh*/
