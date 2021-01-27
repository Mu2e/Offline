#ifndef TrackerConditions_TrackerStatusMaker_hh
#define TrackerConditions_TrackerStatusMaker_hh
//
// Initialize the object defining the complete tracker status
// Initial author: Dave Brown (LBNL) 9/2020
//

// Mu2e includes
#include "TrackerConditions/inc/TrackerStatus.hh"
#include "TrackerConfig/inc/TrackerStatusConfig.hh"
#include "DbTables/inc/TrkElementStatus.hh"

// C++ includes
#include <set>
#include <iosfwd>

namespace mu2e {

  class TrackerStatusMaker {

  public:
    TrackerStatusMaker(TrackerStatusConfig const& config):config_(config) {}
    TrackerStatus::ptr_t fromFcl();
    TrackerStatus::ptr_t fromDb(
      TrkPlaneStatus::cptr_t tpls_p,
      TrkPanelStatus::cptr_t   tpas_p,
      TrkStrawStatusLong::cptr_t   tssl_p,
      TrkStrawStatusShort::cptr_t   tsss_p ); // construct from plane, panel, long-term straw and short-term straw status tables
  private:
    // this object needs to be thread safe, config_ should only be initialized once
    const TrackerStatusConfig config_;
  };

} // namespace mu2e

#endif /* TrackerConditions_TrackerStatusMaker_hh */
