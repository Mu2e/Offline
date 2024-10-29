#ifndef TrackerConditions_TrackerStatus_hh
#define TrackerConditions_TrackerStatus_hh
//
// Define status of the tracker, as the agregate of tracker element status
//

// Mu2e includes
#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/DataProducts/inc/StrawIdMask.hh"
#include "Offline/DataProducts/inc/StrawStatus.hh"
#include "Offline/Mu2eInterfaces/inc/ProditionsEntity.hh"
#include "fhiclcpp/type_traits.h"

// C++ includes
#include <map>
#include <memory>

namespace mu2e {

  class TrackerStatus : public ProditionsEntity {
    public:

      typedef std::shared_ptr<TrackerStatus> ptr_t;
      typedef std::shared_ptr<const TrackerStatus> cptr_t;
      typedef std::map<StrawId,StrawStatus> stat_t;
      typedef std::map<StrawIdMask,stat_t> estat_t;
      constexpr static const char* cxname = {"TrackerStatus"};

      TrackerStatus(): ProditionsEntity(cxname){}

      virtual ~TrackerStatus() {}

      // add status for an element.  This is cumulative
      void addStatus(StrawId const& sid, StrawIdMask const& mask, StrawStatus const& status);

      void print( std::ostream& ) const override;
      // convenience operators for some common situations
      bool noSignal(StrawId const& sid) const;    // return 'true' if we expect no signal from this straw
      bool noisy(StrawId const& sid) const;  // This straw sometimes produces a valid signal, but not always
      bool suppress(StrawId const& sid) const;  // This straw may produce a signal, but it should be suppressed as it is inaccurate
      bool noMaterial(StrawId const& sid) const;  // straw doesn't contribute to scattering or energy loss

      // Net status of an individual Straw.  If the straw is in a plane or panel with status, that will be aggregated
      StrawStatus strawStatus(StrawId const& sid) const;
      // same for panel, plane.  Note; these will return only status that applies to the ENTIRE PANEL (or plane)
      // It will NOT detect (say) a panel where every straw has had the same status individually set (that's not a good configuration)
      StrawStatus panelStatus(StrawId const& sid) const;
      StrawStatus planeStatus(StrawId const& sid) const;
    private:
      std::array<StrawStatus,StrawId::_nustraws> _fullstatus;
      estat_t _status;
  };

} // namespace mu2e

#endif /* TrackerConditions_TrackerStatus_hh */
