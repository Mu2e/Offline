//
// Hold information about the toy tracker, the CTracker.
//
//
// $Id: CTracker.cc,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
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




