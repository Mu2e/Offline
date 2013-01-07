//
// Details common to many straws.
//
// $Id: StrawDetail.cc,v 1.4 2013/01/07 03:59:37 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/01/07 03:59:37 $
//
// Original author Rob Kutschke
//

#include <sstream>

#include "TrackerGeom/inc/StrawDetail.hh"

using namespace std;

namespace mu2e {

  StrawDetail::StrawDetail( int    id,
                            vector<string> const& materialNames,
                            double radius,
                            double thickness,
                            double halfLength,
                            double activeHalfLength,
                            double rwire
                            ):
    _id(id),
    _materialNames(materialNames),
    _radius(radius),
    _thickness(thickness),
    _halfLength(halfLength),
    _activeHalfLength(activeHalfLength),
    _rwire(rwire)
  {
    // Should throw if this fails.
    // if ( _materialNames.size() != 3 ) {}
  }

  // Return G4TUBS parameters outer volume for this straw, includes gas
  // plus straw material.
  TubsParams StrawDetail::getOuterTubsParams() const{
    return TubsParams ( 0., outerRadius(), halfLength() );
  }

  // Return G4TUBS parameters for the straw skin.
  TubsParams StrawDetail::getWallTubsParams() const{
    return TubsParams ( 0., innerRadius(), halfLength() );
  }

  // Return G4TUBS parameters for the wire.
  TubsParams StrawDetail::getWireTubsParams() const{
    return TubsParams ( 0., wireRadius(), halfLength() );
  }

  TubsParams StrawDetail::getGasOnlyTubsParams() const{
    return TubsParams ( wireRadius(), innerRadius(), halfLength() );
  }

  // Construct a string containing the straw Id.
  std::string StrawDetail::name( std::string const& base ) const{
    std::ostringstream os;
    os << base
       << _id;
    return os.str();
  }

  TubsParams StrawDetail::wallMother() const {
    double rIn = _radius - _thickness;
    return TubsParams( rIn, _radius, _halfLength );
  }

  TubsParams StrawDetail::wallOuterMetal() const{
    double rIn = _radius - _outerMetalThickness;
    return TubsParams( rIn, _radius, _halfLength );
  }

  TubsParams StrawDetail::wallCore() const{
    double rIn  = _radius - _thickness + _innerMetal1Thickness + _innerMetal2Thickness;
    double rOut = _radius - _outerMetalThickness;
    return TubsParams( rIn, rOut, _halfLength );
  }

  TubsParams StrawDetail::wallInnerMetal1() const {
    double rIn  = _radius - _thickness + _innerMetal2Thickness;
    double rOut = _radius - _thickness + _innerMetal1Thickness + _innerMetal2Thickness;
    return TubsParams( rIn, rOut, _halfLength );
  }

  TubsParams StrawDetail::wallInnerMetal2() const {
    double rIn  = _radius - _thickness;
    double rOut = _radius - _thickness + _innerMetal2Thickness;
    return TubsParams( rIn, rOut, _halfLength );
  }

  TubsParams StrawDetail::wireMother() const {
    return TubsParams( 0., _rwire, _halfLength );
  }

  TubsParams StrawDetail::wireCore() const {
    double rOut = _rwire-_wirePlateThickness;
    return TubsParams( 0., rOut, _halfLength );
  }

  TubsParams StrawDetail::wirePlate() const {
    double rIn = _rwire-_wirePlateThickness;
    return TubsParams( rIn, _rwire, _halfLength );
  }

} // namespace mu2e
