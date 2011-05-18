#ifndef CTrackerGeom_CTracker_hh
#define CTrackerGeom_CTracker_hh

//
// Hold information about the toy tracker, the CTracker.
// This is just a set of infinitely thin concentric circles.
//
//
// $Id: CTracker.hh,v 1.3 2011/05/18 02:27:14 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:14 $
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
