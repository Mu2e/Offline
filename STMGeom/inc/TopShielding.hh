#ifndef STMGeom_TopShielding_hh
#define STMGeom_TopShielding_hh

// Germanium Detector Object
//
// Author: Haichuan Cao
// Sept 2023

#include <string>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class TopShielding {
  public:
    TopShielding(bool build, bool Skirt_build,
    double T_Zlength,  double T_Xlength, double Front_LT,
    double Al_depth, double Pb_depth, double Cu_depth, double BP_depth,
    double T_ZHole, double BarL, double BarR, double GapL, double GapR, double leak,
    CLHEP::Hep3Vector const & originInMu2e = CLHEP::Hep3Vector(), CLHEP::HepRotation const & rotation = CLHEP::HepRotation()
   ):
      _build(build), _Skirt_build(Skirt_build),
      _T_Zlength(T_Zlength), _T_Xlength(T_Xlength), _Front_LT(Front_LT),
      _Tcontainerdepth(Al_depth), _Tleaddepth(Pb_depth), _Tcopperdepth(Cu_depth), _TBPdepth(BP_depth),
      _T_ZHole(T_ZHole), _Top_Bar_Left(BarL), _Top_Bar_Right(BarR), _Top_Gap_Left(GapL), _Top_Gap_Right(GapR), _Top_Leak(leak),
      _originInMu2e(originInMu2e), _rotation(rotation)
    {
    }

   bool    build()            const {return _build;}
   bool    Skirt_build()      const {return _Skirt_build;}
   double  T_Zlength()        const {return _T_Zlength;}
   double  T_Xlength()        const {return _T_Xlength;}
   double  Front_LT()         const {return _Front_LT;}
   double  Tcontainerdepth()  const {return _Tcontainerdepth;}
   double  Tleaddepth()       const {return _Tleaddepth;}
   double  Tcopperdepth()     const {return _Tcopperdepth;}
   double  TBPdepth()         const {return _TBPdepth;}
   double  T_ZHole()          const {return _T_ZHole;}
   double  Top_Bar_Left()     const {return _Top_Bar_Left;}
   double  Top_Bar_Right()    const {return _Top_Bar_Right;}
   double  Top_Gap_Left()     const {return _Top_Gap_Left;}
   double  Top_Gap_Right()    const {return _Top_Gap_Right;}
   double  Top_Leak()         const {return _Top_Leak;}

   double  Top_Thickness()    const {return _Tleaddepth*2 + _TBPdepth*2 + _Tcopperdepth;}

   TopShielding() {}

  private:

    bool               _build;
    bool               _Skirt_build;
    double             _T_Zlength;
    double             _T_Xlength;
    double             _Front_LT;
    double             _Tcontainerdepth;
    double             _Tleaddepth;
    double             _Tcopperdepth;
    double             _TBPdepth;
    double             _T_ZHole;
    double             _Top_Bar_Left;
    double             _Top_Bar_Right;
    double             _Top_Gap_Left;
    double             _Top_Gap_Right;
    double             _Top_Leak;

    CLHEP::Hep3Vector  _originInMu2e;
    CLHEP::HepRotation _rotation;
  };

}

#endif/*STMGeom_TopShielding_hh*/
