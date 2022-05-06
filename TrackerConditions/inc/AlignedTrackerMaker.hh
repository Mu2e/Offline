#ifndef TrackerConditions_AlignedTrackerMaker_hh
#define TrackerConditions_AlignedTrackerMaker_hh
//
// Make AlignedTracker from fcl or (eventually) database
//

#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/TrackerConfig/inc/AlignedTrackerConfig.hh"
#include "Offline/DbTables/inc/TrkAlignElement.hh"
#include "Offline/DbTables/inc/TrkAlignStraw.hh"

namespace mu2e {


  class AlignedTrackerMaker {
    typedef std::shared_ptr<Tracker> ptr_t;
    public:
    AlignedTrackerMaker(AlignedTrackerConfig const& config):_config(config) {}
    void alignTracker(ptr_t ptr, std::vector<TrkAlignParams> const& tracker_align_params, std::vector<TrkAlignParams> const& plane_align_params, std::vector<TrkAlignParams> const& panel_align_params, std::vector<TrkStrawEndAlign> const& straw_align_params);
    ptr_t fromFcl();
    ptr_t fromDb(TrkAlignTracker::cptr_t tatr_p,
        TrkAlignPlane::cptr_t   tapl_p,
        TrkAlignPanel::cptr_t   tapa_p,
        TrkAlignStraw::cptr_t   tast_p );

    private:

    // this object needs to be thread safe,
    // _config should only be initialized once
    const AlignedTrackerConfig _config;
  };


} // namespace mu2e

#endif /* TrackerConditions_AlignedTrackerMaker_hh */
