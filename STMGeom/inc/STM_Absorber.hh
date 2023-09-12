#ifndef STMGeom_STM_Absorber_hh
#define STMGeom_STM_Absorber_hh

// Germanium Detector Object
//
// Author: Haichuan Cao
// Sept 2023

#include <string>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class STM_Absorber {
  public:
    STM_Absorber(bool build, 
    double Absorber_hW, double Absorber_hH, double Absorber_hT,
    double Absorber_GaptoSSC):
      _build(build),
    _Absorber_hW(Absorber_hW), _Absorber_hH(Absorber_hH), _Absorber_hT(Absorber_hT),
    _Absorber_GaptoSSC(Absorber_GaptoSSC)
    {}

    bool    build()                            const {return _build;}
    double  Absorber_hW()                      const {return _Absorber_hW;}
    double  Absorber_hH()                      const {return _Absorber_hH;}
    double  Absorber_hT()                      const {return _Absorber_hT;}
    double  Absorber_GaptoSSC()                const {return _Absorber_GaptoSSC;}

    STM_Absorber() {}

  private:
    bool               _build;
    double             _Absorber_hW;
    double             _Absorber_hH;
    double             _Absorber_hT;
    double             _Absorber_GaptoSSC;

  };

}

#endif/*STMGeom_STM_Absorber_hh*/
