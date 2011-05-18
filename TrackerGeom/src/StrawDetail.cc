//
// Details common to many straws.
//
// $Id: StrawDetail.cc,v 1.3 2011/05/18 02:27:20 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:20 $
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
                            double rwire
                            ):
    _id(id),
    _materialNames(materialNames),
    _radius(radius),
    _thickness(thickness),
    _halfLength(halfLength),
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


} // namespace mu2e
