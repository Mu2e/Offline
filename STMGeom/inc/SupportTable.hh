#ifndef STMGeom_SupportTable_hh
#define STMGeom_SupportTable_hh

// Table to support things Object
//
// Author: Anthony Palladino
//

#include <string>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class SupportTable {
  public:
    SupportTable(bool build, double tabletopHalfWidth, 
                 double tabletopHalfHeight, double tabletopHalfLength, 
                 double legRadius,
                 CLHEP::Hep3Vector const & originInMu2e = CLHEP::Hep3Vector(),
                 CLHEP::HepRotation const & rotation = CLHEP::HepRotation(), 
                 std::string const & materialName = ""
                ) :
      _build( build ),          
      _tabletopHalfWidth( tabletopHalfWidth ),
      _tabletopHalfHeight( tabletopHalfHeight ),
      _tabletopHalfLength( tabletopHalfLength ),
      _legRadius( legRadius ),
      _originInMu2e( originInMu2e ),
      _rotation    ( rotation     ),
      _materialName( materialName )
    {}

    bool   build()                 const { return _build;  }
    double tabletopHalfWidth()     const { return _tabletopHalfWidth;  }
    double tabletopHalfHeight()     const { return _tabletopHalfHeight;  }
    double tabletopHalfLength()    const { return _tabletopHalfLength; }
    double legRadius()              const { return _legRadius; }
    //double zBegin()          const { return _originInMu2e.z() - zTabletopHalfLength(); }
    //double zEnd()            const { return _originInMu2e.z() + zTabletopHalfLength(); }
   
    CLHEP::Hep3Vector const &  originInMu2e() const { return _originInMu2e; }
    CLHEP::HepRotation const & rotation()     const { return _rotation; }
    std::string const &        materialName() const { return _materialName; }    

    // Genreflex can't do persistency of vector<SupportTable> without a default constructor
    SupportTable() {}
  private:

    bool   _build;
    double _tabletopHalfWidth;
    double _tabletopHalfHeight;
    double _tabletopHalfLength;
    double _legRadius;
    CLHEP::Hep3Vector  _originInMu2e;
    CLHEP::HepRotation _rotation; // wrt to parent volume
    std::string        _materialName;

  };

}

#endif/*STMGeom_SupportTable_hh*/
