#ifndef STMGeom_ElectronicShielding_hh
#define STMGeom_ElectronicShielding_hh

// Germanium Detector Object
//
// Author: Anthony Palladino
//

#include <string>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class ElectronicShielding {
  public:
    ElectronicShielding(bool build, 
    double SiGridX, double SiGridY, double SiGridZ,
    double SiXcenter, double SiYcenter, double SiZcenter,
    double ConcreteT):
      _build(build),
     _SiGridX(SiGridX), _SiGridY(SiGridY), _SiGridZ(SiGridZ),
     _SiXcenter(SiXcenter), _SiYcenter(SiYcenter), _SiZcenter(SiZcenter), 
     _ConcreteT(ConcreteT)    
    {}

    bool    build()                              const {return _build;}
    double  SiGridX()                            const {return _SiGridX;}
    double  SiGridY()                            const {return _SiGridY;}
    double  SiGridZ()                            const {return _SiGridZ;}
    double  SiXcenter()                          const {return _SiXcenter;}
    double  SiYcenter()                          const {return _SiYcenter;}
    double  SiZcenter()                          const {return _SiZcenter;}
    double  ConcreteT()                          const {return _ConcreteT;}

    ElectronicShielding() {}

  private:

    bool               _build;
    double             _SiGridX;
    double             _SiGridY;
    double             _SiGridZ;
    double             _SiXcenter;
    double             _SiYcenter;
    double             _SiZcenter;
    double             _ConcreteT;
   
  };

}

#endif/*STMGeom_ElectronicShielding_hh*/
