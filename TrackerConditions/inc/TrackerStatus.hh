#ifndef TrackerConditions_TrackerStatus_hh
#define TrackerConditions_TrackerStatus_hh
//
// Define status of the tracker, as the agregate of tracker element status
//

// Mu2e includes
#include "DataProducts/inc/StrawId.hh"
#include "DataProducts/inc/StrawIdMask.hh"
#include "TrackerConditions/inc/StrawStatus.hh"
#include "Mu2eInterfaces/inc/ProditionsEntity.hh"
#include "fhiclcpp/type_traits.h"

// C++ includes
#include <set>
#include <memory>

namespace mu2e {

  struct TrackerElementStatus {
    StrawId sid_; // ID of this tracker element.  Only the fields covered by the mask are relevant
    StrawIdMask mask_; // level of this element (straw, panel, or plane)
    StrawStatus status_; // status of this tracker element
    TrackerElementStatus(StrawId const& sid, StrawIdMask const& mask, StrawStatus const& status) : sid_(sid), mask_(mask), status_(status) {}
    bool operator < (const TrackerElementStatus& other) const { return sid_ < other.sid_; } // needed for std::set
  };

  class TrackerStatus : public ProditionsEntity {
  public:

    typedef std::shared_ptr<TrackerStatus> ptr_t;
    typedef std::shared_ptr<const TrackerStatus> cptr_t;
    typedef std::set<TrackerElementStatus> estat_t;
    TrackerStatus(): _name("TrackerStatus") {}
    TrackerStatus(estat_t const& estatus): _name("TrackerStatus"), _estatus(estatus) {}

    virtual ~TrackerStatus() {}

    std::string const& name() const override { return _name; }
    void print( std::ostream& ) const override;
    // convenience operators for some common situations
    bool noSignal(StrawId const& sid) const;    // return 'true' if we expect no usable signal from this straw
    bool suppress(StrawId const& sid) const;  // This straw may produce a signal, but it should be suppressed as it is inaccurate 
    bool noMaterial(StrawId const& sid) const;  // starw doesn't contribut to scattering or energy loss

    // Net status of an individual Straw.  If the straw is in a plane or panel with status, that will be aggregated
    StrawStatus strawStatus(StrawId const& sid) const;
    // same for panel, plane.  Note; these will return only status that applies to the ENTIRE PANEL (or plane)
    // It will NOT detect (say) a panel where every straw has had the same status individually set (that's not a good configuration)
    StrawStatus panelStatus(StrawId const& sid) const;
    StrawStatus planeStatus(StrawId const& sid) const;
  private:
    std::string _name;
    estat_t _estatus; // sparse list of tracker element status
  };

} // namespace mu2e

#endif /* TrackerConditions_TrackerStatus_hh */
