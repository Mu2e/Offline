#ifndef ExampleExtrasCTracker_HH
#define ExampleExtrasCTracker_HH

//
// Hold information about the toy tracker, the CTracker.
// This is just a set of infinitely thin concentric circles.
//
//
// $Id: CTracker.hh,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
//
// Original author Rob Kutschke
//

#include <vector>
#include "GeometryService/inc/Detector.hh"

namespace mu2e {

// Forward reference.
class SimpleConfig;

  class CTracker: public Detector{
    
 public:
  CTracker( const SimpleConfig& c );
  ~CTracker();
  
  // Accessors
  const std::vector<double>& r() const { return radii;}
  double zHalf()const { return zhalf;}
    
 private:
  std::vector<double> radii;
  double zhalf;
  
  };

}

#endif
