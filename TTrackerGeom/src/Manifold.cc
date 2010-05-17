//#include <sstream>

#include "TTrackerGeom/inc/Manifold.hh"

using namespace std;

using CLHEP::Hep3Vector;

namespace mu2e {

  // Default constructor
  Manifold::Manifold():
    _origin(CLHEP::Hep3Vector(0.,0.,0.)){
    _halfLengths.push_back(-1.);
    _halfLengths.push_back(-1.);
    _halfLengths.push_back(-1.);
  }

  // Main constructor
  Manifold::Manifold(const CLHEP::Hep3Vector& origin,
                     const vector<double>& halfLengths
                     ):
    _origin(origin),
    _halfLengths(halfLengths)
  {}
  
  // Construct a string containing the device Id.
  /*
    std::string Manifold::name( std::string const& base ) const{
    std::ostringstream os;
    os << base 
    << _id.getSectorId()._did    << "_"
         << _id.getSectorId()._sector << "_"
         << _id.getManifold();
         return os.str();
         }
  */
  
} // end namespace mu2e
