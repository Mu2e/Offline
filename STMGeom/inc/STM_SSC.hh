#ifndef STMGeom_STM_SSC_hh
#define STMGeom_STM_SSC_hh

// STM Collimator Object
//
// Author: Haichuan Cao
// Sept 2023

#include <string>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class STM_SSC {
  public:

    STM_SSC(bool build, bool VDbuild,
                 double delta_WlR, double delta_WlL, double W_middle,
                 double W_height, double Wdepth_f, double Wdepth_b,
                 double Aperture_HPGe1, double Aperture_HPGe2, double Aperture_LaBr1, double Aperture_LaBr2,
                 double offset_Spot, double leak, double FrontToWall, double ZGap, double ZGapBack,
                 CLHEP::Hep3Vector const & originInMu2e = CLHEP::Hep3Vector(),
                 CLHEP::HepRotation const & rotation = CLHEP::HepRotation(),
                 std::string const & material = ""
                 ) :
      _build(build),
      _VDbuild(build),
      _delta_WlR(delta_WlR),
      _delta_WlL(delta_WlL),
      _W_middle(W_middle),
      _W_height(W_height),
      _Wdepth_f(Wdepth_f),
      _Wdepth_b(Wdepth_b),
      _Aperture_HPGe1(Aperture_HPGe1),
      _Aperture_HPGe2(Aperture_HPGe2),
      _Aperture_LaBr1(Aperture_LaBr1),
      _Aperture_LaBr2(Aperture_LaBr2),
      _offset_Spot(offset_Spot),
      _leak(leak),
      _FrontToWall(FrontToWall),
      _ZGap(ZGap),
      _ZGapBack(ZGapBack),
      _originInMu2e(originInMu2e),
      _rotation(rotation),
      _material(material)
    {
    }

    bool   build()       const {return _build;}
    bool   VDbuild()     const {return _VDbuild;}
    double delta_WlR()   const {return _delta_WlR;}
    double delta_WlL()   const {return _delta_WlL;}
    double W_middle()    const {return _W_middle;}
    double W_length()    const {return _W_middle+_delta_WlR+_delta_WlL;}
    double W_height()    const {return _W_height;}
    double Wdepth_f()    const {return _Wdepth_f;}
    double Wdepth_b()    const {return _Wdepth_b;}
    double W_depth()     const {return _Wdepth_f+_Wdepth_b;}

    double Aperture_HPGe1()    const {return _Aperture_HPGe1;}
    double Aperture_HPGe2()    const {return _Aperture_HPGe2;}
    double Aperture_LaBr1()    const {return _Aperture_LaBr1;}
    double Aperture_LaBr2()    const {return _Aperture_LaBr2;}
    double offset_Spot()       const {return _offset_Spot;}
    double leak()              const {return _leak;}
    double FrontToWall()       const {return _FrontToWall;}
    double ZGap()              const {return _ZGap;}
    double ZGapBack()          const {return _ZGapBack;}


    //double zBegin()          const { return _originInMu2e.z() - zTabletopHalfLength(); }
    //double zEnd()            const { return _originInMu2e.z() + zTabletopHalfLength(); }

    CLHEP::Hep3Vector const &  originInMu2e()     const { return _originInMu2e; }
    CLHEP::HepRotation const & rotation()         const { return _rotation; }
    std::string const &        material()         const { return _material; }
    // Genreflex can't do persistency of vector<STM_SSC> without a default constructor
    STM_SSC() {}

  private:

    bool   _build;
    bool   _VDbuild;
    double _delta_WlR;
    double _delta_WlL;
    double _W_middle;
    double _W_height;
    double _Wdepth_f;
    double _Wdepth_b;
    double _Aperture_HPGe1;
    double _Aperture_HPGe2;
    double _Aperture_LaBr1;
    double _Aperture_LaBr2;
    double _offset_Spot;
    double _leak;
    double _FrontToWall;
    double _ZGap;
    double _ZGapBack;

    CLHEP::Hep3Vector  _originInMu2e;
    CLHEP::HepRotation _rotation; // wrt to parent volume
    std::string        _material;
  };

}

#endif/*STMGeom_STM_SSC_hh*/
