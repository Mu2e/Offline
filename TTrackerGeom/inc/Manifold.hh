#ifndef TTRACKER_MANIFOLD_HH
#define TTRACKER_MANIFOLD_HH

#include <vector>
#include <string>

#include "TTrackerGeom/inc/ManifoldId.hh"

#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class Manifold{

  public:
    
    // Default constructor
    Manifold();
    
    // Main constructor
    Manifold( const CLHEP::Hep3Vector& origin,
              const std::vector<double>& halfLengths
              );
    
    // Destructor
    ~Manifold(){};
    
    // Return origin in tracker coordinates
    const CLHEP::Hep3Vector& getOrigin() const {return _origin;}

    //Return origin in coordinates local to the device
    //const CLHEP::Hep3Vector& getDeviceLocalOrigin() const {
    //return _origin;
    //}
    
    // Return halfLengths
    const std::vector<double>& getHalfLengths() const {return _halfLengths;}
    
    // Construct a string containing the device Id.
    //std::string name( std::string const& base ) const;

  protected:
    
    // Member Variables
    
    CLHEP::Hep3Vector   _origin;
    std::vector<double> _halfLengths;
    
  };
  
} // end namespace mu2e

#endif
