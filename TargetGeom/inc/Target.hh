#ifndef ExampleExtrasTarget_HH
#define ExampleExtrasTarget_HH

//
// Class to represent the system of target foils.
// For now these are just disks perpendicular to the z axis.
//
// $Id: Target.hh,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
//
// Original author Rob Kutschke
//
// Coordinates are given in the detector coordinate 
// system in cm.
//

// Includes from C++
#include <vector>

// Includes from Mu2e
#include "TargetGeom/inc/TargetFoil.hh"
#include "GeometryService/inc/Detector.hh"

namespace mu2e {

  // Forward reference.
  class SimpleConfig;

  class Target : public Detector{

  public:
    Target( const SimpleConfig& c);
    ~Target();
    
    int nFoils() const { return _foils.size(); }

    TargetFoil const& foil( int n ) const { return _foils.at(n); }
    
  private:

    // All dimensions in cm.
    std::vector<TargetFoil> _foils;
    
};

}
#endif
