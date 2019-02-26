#ifndef TrackerGeom_StrawDetail_hh
#define TrackerGeom_StrawDetail_hh

//
// Class to hold information about the properties of each type of straw.
// We need different types of straws for the different lengths of
// straws in the Tracker.
//
//
// $Id: StrawDetail.hh,v 1.10 2013/03/26 23:28:23 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/26 23:28:23 $
//
// Original author Rob Kutschke
//

#include "GeomPrimitives/inc/TubsParams.hh"

#include <string>
#include <vector>

namespace mu2e {

  class StrawDetail{

    friend class TrackerMaker;

  public:
    StrawDetail():
      _id(-1),
      _materialNames(),
      _radius(-1.),
      _thickness(-1.),
      _halfLength(-1.),
      _activeHalfLength(-1.),
      _rwire(-1.)
    {}

    StrawDetail( int    id,
                 std::vector<std::string> const& materialNames,
                 double radius,
                 double thickness,
                 double halfLength,
                 double activehalfLength,
                 double rwire
                 );

    // Compiler generated versions are OK for destructor
    // and for copy and assignment constructors.

    int Id() const { return  _id; }

    double outerRadius()     const { return _radius;}
    double innerRadius()     const { return _radius-_thickness;}
    double wireRadius()      const { return _rwire; }
    double thickness()       const { return _thickness; }
    double length()          const { return _halfLength*2.0; }
    double halfLength()      const { return _halfLength; }
    double activeHalfLength() const { return _activeHalfLength; }
    double activeLength()     const { return _activeHalfLength*2.0; }

    std::string const& wallMaterialName() const{ return _materialNames[0]; }
    std::string const&  gasMaterialName() const{ return _materialNames[1]; }
    std::string const& wireMaterialName() const{ return _materialNames[2]; }


    // Return G4TUBS parameters outer volume for this straw, includes
    // wire, gas and straw materials.
    TubsParams getOuterTubsParams() const;

    // Return G4TUBS parameters for the straw skin.
    TubsParams getWallTubsParams() const;

    // Return G4TUBS parameters for the wire.
    TubsParams getWireTubsParams() const;

    // Return G4TUBS parameters for the gas volume proper.
    TubsParams getGasOnlyTubsParams() const;

    std::string const& materialName(int idx) const { return _materialNames.at(idx);}

    std::vector<std::string> materialNames() const { return _materialNames;}

    std::string name( std::string const& base ) const;

    // The next two groups are defined only for SupportModel::detailedv0
    TubsParams wallMother()      const;
    TubsParams wallOuterMetal()  const;
    TubsParams wallCore()        const;
    TubsParams wallInnerMetal1() const;
    TubsParams wallInnerMetal2() const;
    TubsParams wireMother()      const;
    TubsParams wireCore()        const;
    TubsParams wirePlate()       const;

    std::string const& wallMotherMaterialName()      const{ return  wallMaterialName();  }
    std::string const& wallOuterMetalMaterialName()  const{ return _outerMetalMaterial;  }
    std::string const& wallCoreMaterialName()        const{ return  wallMaterialName();  }
    std::string const& wallInnerMetal1MaterialName() const{ return _innerMetal1Material; }
    std::string const& wallInnerMetal2MaterialName() const{ return _innerMetal2Material; }
    std::string const& wireMotherMaterialName()      const{ return  wireMaterialName();  }
    std::string const& wireCoreMaterialName()        const{ return  wireMaterialName();  }
    std::string const& wirePlateMaterialName()       const{ return _wirePlateMaterial;   }


  private:

    // Identifier for this type of straw.
    int _id;

    // Order of materials is:
    // straw material, gas volume, wire
    std::vector<std::string> _materialNames;

    double _radius;            // Outer radius of straw, gas
    double _thickness;         // Thickness of the straw material
    double _halfLength;        // Physical length of the straw
    double _activeHalfLength;  // This may be shorter than the physical length if the ends are deadened.
    double _rwire;             // Radius of the wire, including plating

    // These are defined only for SupportModel::detailedv0
    double _outerMetalThickness;  // Outer metal coating on the straw
    double _innerMetal1Thickness; // Metal coating immediately on the inside of the straw.
    double _innerMetal2Thickness; // Metal coating on the inside of Metal1
    double _wirePlateThickness;   // Plating on the wire.
    std::string _outerMetalMaterial;  // Outer metal coating on the straw
    std::string _innerMetal1Material; // Metal coating immediately on the inside of the straw.
    std::string _innerMetal2Material; // Metal coating on the inside of Metal1
    std::string _wirePlateMaterial;   // Plating on the wire.

  };

}  //namespace mu2e

#endif /* TrackerGeom_StrawDetail_hh */
