#ifndef TrackerConditions_AlignedTrackerMaker_hh
#define TrackerConditions_AlignedTrackerMaker_hh
//
// Make AlignedTracker from fcl or (eventually) database
//

#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerConditions/inc/AlignedTrackerConfig.hh"
#include "DbTables/inc/TrkAlignTracker.hh"
#include "DbTables/inc/TrkAlignPlane.hh"
#include "DbTables/inc/TrkAlignPanel.hh"

namespace mu2e {

  class AlignedTrackerMaker {

  public:
    AlignedTrackerMaker(AlignedTrackerConfig const& config):_config(config) {}
    Tracker::ptr_t fromFcl();
    Tracker::ptr_t fromDb(TrkAlignTracker::cptr_t tatr_p,
			  TrkAlignPlane::cptr_t   tapl_p,
			  TrkAlignPanel::cptr_t   tapa_p );

  private:

    // this object needs to be thread safe, 
    // _config should only be initialized once
    const AlignedTrackerConfig _config;
  };


} // namespace mu2e

#endif /* TrackerConditions_AlignedTrackerMaker_hh */
