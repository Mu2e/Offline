#ifndef TrackerGeom_StrawDetail_hh
#define TrackerGeom_StrawDetail_hh

//
// Class to hold information about the properties of each type of straw.
// We need different types of straws for the conducting and non-conducting
// straws in the LTracker and for different lengths of straws in the TTracker.
//
//
// $Id: StrawDetail.hh,v 1.5 2011/05/18 02:27:19 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:19 $
//
// Original author Rob Kutschke
//

#include "TrackerGeom/inc/TubsParams.hh"

#include <vector>
#include <string>

namespace mu2e {

  class StrawDetail{

  public:
    StrawDetail():
      _id(-1),
      _materialNames(),
      _radius(-1.),
      _thickness(-1.),
      _halfLength(-1.),
      _rwire(-1.)
    {}

    StrawDetail( int    id,
                 std::vector<std::string> const& materialNames,
                 double radius,
                 double thickness,
                 double length,
                 double rwire
                 );

    // Compiler generated versions are OK for destructor
    // and for copy and assignment constructors.

    int Id() const { return  _id; }

    double outerRadius()   const { return _radius;}
    double innerRadius()   const { return _radius-_thickness;}
    double wireRadius()    const { return _rwire; }
    double thickness()     const { return _thickness; }
    double length()        const { return _halfLength*2.0; }
    double halfLength()    const { return _halfLength; }

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

  private:

    // Identifier for this type of straw.
    int32_t _id;

    // Order of materials is:
    // straw material, gas volume, wire
    std::vector<std::string> _materialNames;

    double _radius;
    double _thickness;
    double _halfLength;
    double _rwire;

  };

}  //namespace mu2e

#endif /* TrackerGeom_StrawDetail_hh */
