#ifndef CTrackerGeom_CTracker_hh
#define CTrackerGeom_CTracker_hh

//
// Hold information about the toy tracker, the CTracker.
// This is just a set of infinitely thin concentric circles.
//
//
// $Id: CTracker.hh,v 1.2 2011/05/17 15:41:35 greenc Exp $
// $Author: greenc $ 
// $Date: 2011/05/17 15:41:35 $
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

#endif /* CTrackerGeom_CTracker_hh */
