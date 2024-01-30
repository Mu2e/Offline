//
// Alignment tables for tracker elements for simulation
// derived from the reco version of the tables, and have the same structure
//

#ifndef DbTables_TrkAlignElementSim_hh
#define DbTables_TrkAlignElementSim_hh

#include "Offline/DbTables/inc/TrkAlignElement.hh"

namespace mu2e {


// unique classes for Db usage: not sure why this is needed
class TrkAlignPanelSim : public TrkAlignPanel {
 public:
  constexpr static const char* cxname = "TrkAlignPanelSim";
  TrkAlignPanelSim() :
      TrkAlignPanel(cxname, "trk.alignpanelsim") {}
};

class TrkAlignPlaneSim : public TrkAlignPlane {
 public:
  constexpr static const char* cxname = "TrkAlignPlaneSim";
  TrkAlignPlaneSim() :
      TrkAlignPlane(cxname, "trk.alignplanesim") {}
};

class TrkAlignTrackerSim : public TrkAlignTracker {
 public:
  constexpr static const char* cxname = "TrkAlignTrackerSim";
  TrkAlignTrackerSim() : TrkAlignTracker(cxname, "trk.aligntrackersim") {}
};

}  // namespace mu2e
#endif
