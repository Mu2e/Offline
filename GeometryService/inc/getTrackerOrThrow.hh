#ifndef GeometryService_getTrackerOrThrow_hh
#define GeometryService_getTrackerOrThrow_hh
//
// Look in the geometry service to find a tracker that
// derives from the Tracker base class ( LTracker or
// TTracker ). Return a reference to the one that is present.
// If neither are present, throw.
//
// $Id: getTrackerOrThrow.hh,v 1.2 2011/05/17 15:41:35 greenc Exp $
// $Author: greenc $ 
// $Date: 2011/05/17 15:41:35 $
//
// Original author Rob Kutschke
//


namespace mu2e {
  class Tracker;

  const Tracker& getTrackerOrThrow();

}
#endif /* GeometryService_getTrackerOrThrow_hh */
