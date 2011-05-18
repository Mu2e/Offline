//
// Hold information about the toy tracker, the CTracker.
//
//
// $Id: CTracker.cc,v 1.2 2011/05/18 02:27:15 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:15 $
//
// Original author Rob Kutschke
//

#include "CTrackerGeom/inc/CTracker.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"

namespace mu2e{

CTracker::CTracker( const SimpleConfig& c ){
  c.getVectorDouble("ctracker.radii",radii);
  zhalf = c.getDouble("ctracker.halfZlength");
}

CTracker::~CTracker(){}

}




