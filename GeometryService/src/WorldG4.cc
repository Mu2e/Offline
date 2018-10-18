// G4 specific geometry info.  Not available in non-geant jobs.
//
// Rob Kutschke, 2011

#include "GeometryService/inc/WorldG4.hh"

#include <limits>
#include <algorithm>
#include <iterator>

#include "cetlib_except/exception.h"

using namespace std;

namespace mu2e {

  WorldG4::WorldG4()
    : _dirtG4Ymin(std::numeric_limits<double>::max())
    , _dirtG4Ymax(std::numeric_limits<double>::min())
  {}

  // The argument is a point in Mu2e coordinates.  Return true if this point
  // is inside the G4 world, false otherwise.
  bool WorldG4::inWorld( CLHEP::Hep3Vector const& x0InMu2e ) const{

    CLHEP::Hep3Vector x0InG4 = x0InMu2e + _mu2eOriginInWorld;
    if ( std::abs(x0InG4.x()) > _halfLengths[0] ) return false;
    if ( std::abs(x0InG4.y()) > _halfLengths[1] ) return false;
    if ( std::abs(x0InG4.z()) > _halfLengths[2] ) return false;

    return true;
  }

  // The arugment is a point in Mu2e coordinates.  Throw an exception if this point
  // is outside the G4 world.
  void WorldG4::inWorldOrThrow( CLHEP::Hep3Vector const& x0 ) const{

    if ( ! inWorld(x0) ){
      throw cet::exception("GEOM")
        << "Point is outside of the G4 world: " << x0 << "\n";
    }
    return;

  }

  std::ostream& operator<<(std::ostream& os, const WorldG4& w) {
    os<<"WordG4(halfLengths = {";
    std::copy(w.halfLengths().begin(), w.halfLengths().end(), std::ostream_iterator<double>(os, ", "));
    os<<"}, hallFormalHalfSize = {";
    std::copy(w.hallFormalHalfSize().begin(), w.hallFormalHalfSize().end(), std::ostream_iterator<double>(os, ", "));
    os<<"}, hallFormalCenterInWorld = "<<w.hallFormalCenterInWorld()
      <<", mu2eOriginInWorld = "<<w.mu2eOriginInWorld()
      <<")";

    return os;
  }

}
