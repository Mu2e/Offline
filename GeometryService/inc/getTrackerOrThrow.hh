#ifndef GeometryService_getTrackerOrThrow_hh
#define GeometryService_getTrackerOrThrow_hh
//
// Look in the geometry service to find a tracker that
// derives from the Tracker base class ( LTracker or
// TTracker ). Return a reference to the one that is present.
// If neither are present, throw.
//
// $Id: getTrackerOrThrow.hh,v 1.3 2011/05/18 02:27:16 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:16 $
//
// Original author Rob Kutschke
//


namespace mu2e {
  class Tracker;

  const Tracker& getTrackerOrThrow();

}
#endif /* GeometryService_getTrackerOrThrow_hh */
