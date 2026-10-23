#ifndef STMGeom_STM_SSC_hh
#define STMGeom_STM_SSC_hh

// STM Collimator Object
//
// Author: Haichuan Cao
// Sept 2023
// Update: Yongyi Wu Sept-Oct 2026, updated (hand-stacked) SSC

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
      // Was _VDbuild(build), which ignored VDbuild. Every geometry file
      // sets the two alike (both true), so no file's result changes.
      _VDbuild(VDbuild),
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
      _material(material),
      _r_LaBr_f(0.),
      _r_HPGe_f(0.),
      _r_LaBr_b(0.),
      _r_HPGe_b(0.),
      _boreToBase(0.)
    {
    }

    // Updated (hand-stacked) SSC.
    //
    // The updated collimator is a single symmetric tungsten block, 6.4 x 3.2
    // x 4.0 in, rather than a middle block flanked by two wings of unequal
    // width. So there is no W_middle/delta_WlL/delta_WlR to give: the total
    // width IS W_width, and the wing-asymmetry offset (delta_WlL-delta_WlR)/2
    // vanishes.
    //
    // The bores are stepped in z, at x = +/-offset_Spot, y = 0. They are given
    // as radii rather than as apertures.
    //
    // Taking a separate constructor rather than defaulted arguments keeps
    // each description to the parameters it actually has, so neither can be
    // built with a quantity that does not apply to it.
    STM_SSC(bool build, bool VDbuild,
                 double W_width, double W_height, double Wdepth_f, double Wdepth_b,
                 double r_LaBr_f, double r_HPGe_f, double r_LaBr_b, double r_HPGe_b,
                 double offset_Spot, double leak, double FrontToWall, double ZGap, double ZGapBack,
                 double boreToBase,
                 CLHEP::Hep3Vector const & originInMu2e = CLHEP::Hep3Vector(),
                 CLHEP::HepRotation const & rotation = CLHEP::HepRotation(),
                 std::string const & material = ""
                 ) :
      _build(build),
      _VDbuild(VDbuild),
      _delta_WlR(0.),
      _delta_WlL(0.),
      _W_middle(W_width),
      _W_height(W_height),
      _Wdepth_f(Wdepth_f),
      _Wdepth_b(Wdepth_b),
      _Aperture_HPGe1(0.),
      _Aperture_HPGe2(0.),
      _Aperture_LaBr1(0.),
      _Aperture_LaBr2(0.),
      _offset_Spot(offset_Spot),
      _leak(leak),
      _FrontToWall(FrontToWall),
      _ZGap(ZGap),
      _ZGapBack(ZGapBack),
      _originInMu2e(originInMu2e),
      _rotation(rotation),
      _material(material),
      _r_LaBr_f(r_LaBr_f),
      _r_HPGe_f(r_HPGe_f),
      _r_LaBr_b(r_LaBr_b),
      _r_HPGe_b(r_HPGe_b),
      _boreToBase(boreToBase)
    {
      // W_width is stored as _W_middle with both wings zero, so the inherited
      // W_length() below returns the full width for either description and
      // the leak VDs and front shielding need no branch of their own.
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

    // Updated (hand-stacked) SSC: the total width, and the stepped bore
    // radii. Zero for the earlier description, which uses W_middle plus the
    // two wings and the Aperture_* areas above instead.
    double W_width()           const {return _W_middle;}
    double r_LaBr_f()          const {return _r_LaBr_f;}
    double r_HPGe_f()          const {return _r_HPGe_f;}
    double r_LaBr_b()          const {return _r_LaBr_b;}
    double r_HPGe_b()          const {return _r_HPGe_b;}

    // Height of the bore axis above the baseplate. Every hand-stacked
    // piece downstream measures its own height up from that baseplate,
    // so this is the datum for the whole shield house.
    double boreToBase()        const {return _boreToBase;}


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

    // updated (hand-stacked) SSC only; zero for the earlier description
    double _r_LaBr_f;
    double _r_HPGe_f;
    double _r_LaBr_b;
    double _r_HPGe_b;
    double _boreToBase;
  };

}

#endif/*STMGeom_STM_SSC_hh*/
