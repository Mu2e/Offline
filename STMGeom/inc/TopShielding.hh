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
    double LiftBeam_L, double LiftBeam_H, double LiftBeam_T, double LiftBeam_Xmove,
    double T_Zlength, double T_Xlength, double Front_LT,
    double TTF_Zlength, double TTF_Xlength, double TTB_Zlength,
    double Al_depth, double Pb_depth, double Cu_depth, double BP_depth,
    double T_ZHole, double BarL, double BarR, double GapL, double GapR, double leak,
    CLHEP::Hep3Vector const & originInMu2e = CLHEP::Hep3Vector(), CLHEP::HepRotation const & rotation = CLHEP::HepRotation()
   ):
      _build(build), _Skirt_build(Skirt_build),
      _LiftBeam_L(LiftBeam_L), _LiftBeam_H(LiftBeam_H), _LiftBeam_T(LiftBeam_T), _LiftBeam_Xmove(LiftBeam_Xmove),
      _T_Zlength(T_Zlength), _T_Xlength(T_Xlength), _Front_LT(Front_LT),
      _TTF_Zlength(TTF_Zlength), _TTF_Xlength(TTF_Xlength), _TTB_Zlength(TTB_Zlength),
      _Tcontainerdepth(Al_depth), _Tleaddepth(Pb_depth), _Tcopperdepth(Cu_depth), _TBPdepth(BP_depth),
      _T_ZHole(T_ZHole), _Top_Bar_Left(BarL), _Top_Bar_Right(BarR), _Top_Gap_Left(GapL), _Top_Gap_Right(GapR), _Top_Leak(leak),
      _originInMu2e(originInMu2e), _rotation(rotation)
    {
    }

   bool    build()            const {return _build;}
   bool    Skirt_build()      const {return _Skirt_build;}

   double  LiftBeam_L()        const {return _LiftBeam_L;}
   double  LiftBeam_H()        const {return _LiftBeam_H;}
   double  LiftBeam_T()        const {return _LiftBeam_T;}
   double  LiftBeam_Xmove()    const {return _LiftBeam_Xmove;}


   double  T_Zlength()        const {return _T_Zlength;}
   double  T_Xlength()        const {return _T_Xlength;}
   double  Front_LT()         const {return _Front_LT;}

   double  TTF_Zlength()      const {return _TTF_Zlength;}
   double  TTF_Xlength()      const {return _TTF_Xlength;}
   double  TTB_Zlength()      const {return _TTB_Zlength;}

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

    double             _LiftBeam_L;
    double             _LiftBeam_H;
    double             _LiftBeam_T;
    double             _LiftBeam_Xmove;

    double             _T_Zlength;
    double             _T_Xlength;
    double             _Front_LT;

    double             _TTF_Zlength;
    double             _TTF_Xlength;
    double             _TTB_Zlength;

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
