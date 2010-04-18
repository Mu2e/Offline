#ifndef getTrackerOrThrow_HH
#define getTrackerOrThrow_HH
//
// Look in the geometry service to find a tracker that
// derives from the Tracker base class ( LTracker or
// TTracker ). Return a reference to the one that is present.
// If neither are present, throw.
//
// $Id: getTrackerOrThrow.hh,v 1.1 2010/04/18 00:32:42 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/04/18 00:32:42 $
//
// Original author Rob Kutschke
//


namespace mu2e {
  class Tracker;

  const Tracker& getTrackerOrThrow();

}
#endif
