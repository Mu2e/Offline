#ifndef MECOStyleProtonAbsorberGeom_MECOStyleProtonAbsorberPart_hh
#define MECOStyleProtonAbsorberGeom_MECOStyleProtonAbsorberPart_hh

//
// Class to represent one part of MECO style Proton Absorber
//
//
// Original author MyeongJae Lee
//
//

#include <string>

// Includes from CLHEP
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class MECOStyleProtonAbsorberPart {

  public:
    MECOStyleProtonAbsorberPart( int id,
                                 CLHEP::Hep3Vector const& c,
                                 double rOut0,
                                 double rIn0,
                                 double rOut1,
                                 double rIn1,
                                 double halflen,
                                 std::string const& m):
      _id(id),
      _c(c),
      _rOut0(rOut0),
      _rIn0(rIn0),
      _rOut1(rOut1),
      _rIn1(rIn1),
      _halflen(halflen),
      _nSides(0),
      _material(m){
    }

    MECOStyleProtonAbsorberPart( int id,
                                 CLHEP::Hep3Vector const& c,
                                 double rOut0,
                                 double rIn0,
                                 double rOut1,
                                 double rIn1,
                                 double halflen,
                                 int    nSides,
                                 std::string const& m):
      _id(id),
      _c(c),
      _rOut0(rOut0),
      _rIn0(rIn0),
      _rOut1(rOut1),
      _rIn1(rIn1),
      _halflen(halflen),
      _nSides(nSides),
      _material(m){
    }

    // Use compiler-generated copy c'tor, copy assignment, and d'tor

    int id() const { return _id; }

    CLHEP::Hep3Vector const& center()  const { return _c;}

    double outerRadiusAtStart()        const { return _rOut0;}
    double innerRadiusAtStart()        const { return _rIn0;}
    double outerRadiusAtEnd()          const { return _rOut1;}
    double innerRadiusAtEnd()          const { return _rIn1;}
    double halfLength()                const { return _halflen;}

    int    nSides()                    const { return _nSides; }
    std::string material()             const { return _material;}


  private:

    // 0 for DS2 part, 1 for DS3 part
    int _id;

    // Center of the part in mu2e coordinate
    CLHEP::Hep3Vector _c;

    // Inner and outer radii at start
    double _rOut0;
    double _rIn0;

    // Inner and outer radii at end
    double _rOut1;
    double _rIn1;

    // halflength
    double _halflen;

    int _nSides;  // 0 for cones, non-zero for barrels with slats

    std::string _material;


  };

}
#endif /* MECOStyleProtonAbsorberGeom_MECOStyleProtonAbsorberPart_hh */
