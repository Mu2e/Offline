#ifndef STMGeom_FrontShielding_hh
#define STMGeom_FrontShielding_hh

// Germanium Detector Object
//
// Author: Haichuan Cao
// Sept 2023

#include <string>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class FrontShielding {
  public:
    FrontShielding(bool build, double height, double W_length, double W_depth,
    double Pb_depth1, double Pb_depth2, double Al_depth, double Cu_depth, double BP_depth,
    double fPb_lengthL, double fPb_lengthR, double GapForTop, double LeakForSSC, double CopperL, double Front_H, double FrontHole_r,
    CLHEP::Hep3Vector const & originInMu2e = CLHEP::Hep3Vector(), CLHEP::HepRotation const & rotation = CLHEP::HepRotation()
    ):
      _build(build), _HeightofRoom(height), _Ftungstenlength(W_length), _Ftungstendepth(W_depth),
      _Fleaddepth1(Pb_depth1), _Fleaddepth2(Pb_depth2), _Faluminumdepth(Al_depth), _Fcopperdepth(Cu_depth), _FBPdepth(BP_depth),
      _fPb_lengthL(fPb_lengthL), _fPb_lengthR(fPb_lengthR), _GapForTop(GapForTop), _LeakForSSC(LeakForSSC), _CopperL(CopperL), _Front_H(Front_H), _FrontHole_r(FrontHole_r),
      _originInMu2e(originInMu2e), _rotation(rotation)
    {}

   bool  build()                     const { return _build;}

   double  HeightofRoom()     const {return _HeightofRoom;}
   double  Ftungstenlength()  const {return _Ftungstenlength;}
   double  Ftungstendepth()   const {return _Ftungstendepth;}
   double  Fleaddepth1()      const {return _Fleaddepth1;}
   double  Fleaddepth2()      const {return _Fleaddepth2;}
   double  Faluminumdepth()   const {return _Faluminumdepth;}
   double  Fcopperdepth()     const {return _Fcopperdepth;}
   double  FBPdepth()         const {return _FBPdepth;}
   double  fPb_lengthL()      const {return _fPb_lengthL;}
   double  fPb_lengthR()      const {return _fPb_lengthR;}
   double  GapForTop()        const {return _GapForTop;}
   double  LeakForSSC()       const {return _LeakForSSC;}
   double  CopperL()          const {return _CopperL;}
   double  Front_H()          const {return _Front_H;}
   double  FrontHole_r()      const {return _FrontHole_r;}

   double  Front_Thickness()  const {return _Ftungstendepth + _Fleaddepth2*3 + _FBPdepth*2;}
   double  Front_Length()     const {return _Ftungstenlength + 2*_LeakForSSC + _fPb_lengthL + _fPb_lengthR;}

    FrontShielding() {}
  private:

    bool               _build;
    double             _HeightofRoom;
    double             _Ftungstenlength;
    double             _Ftungstendepth;
    double             _Fleaddepth1;
    double             _Fleaddepth2;
    double             _Fleaddepth4;
    double             _Faluminumdepth;
    double             _Fcopperdepth;
    double             _FBPdepth;
    double             _fPb_lengthL;
    double             _fPb_lengthR;
    double             _GapForTop;
    double             _LeakForSSC;
    double             _CopperL;
    double             _Front_H;
    double             _FrontHole_r;

    CLHEP::Hep3Vector  _originInMu2e;
    CLHEP::HepRotation _rotation;
  };

}

#endif/*STMGeom_FrontShielding_hh*/
