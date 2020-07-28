#ifndef TrackerConditions_TrackerStatus_hh
#define TrackerConditions_TrackerStatus_hh
//
// Define status of the tracker
//

// Mu2e includes
#include "DataProducts/inc/StrawId.hh"
#include "DataProducts/inc/StrawIdMask.hh"
#include "TrackerConditions/inc/StrawStatus.hh"
#include "Mu2eInterfaces/inc/ProditionsEntity.hh"

// C++ includes
#include <set>
#include <memory>

namespace mu2e {

// struct to define a range of deadness in a straw.  This is assumed to be equidistant
// from the straw center, to model radiation damage due to dose (which is concentrated
// at the lowest radius = straw center)
//

  class TrackerStatus : public ProditionsEntity {
  public:
    Struct TrackerElementStatus {
      StrawId sid_; // ID of this tracker element
      StrawIdMask mask_; // level of this element (straw, panel, or plane)
      StrawStatus status_; // status of this tracker element
    };

    typedef std::shared_ptr<TrackerStatus> ptr_t;
    typedef std::shared_ptr<const TrackerStatus> cptr_t;
    typedef std::set<TrackerElementStatus> estat_t;
    typedef std::set<StrawStatus> sstat_t;
    TrackerStatus(): _name("TrackerStatus") {}
    TrackerStatus(set_t const& deadstraws):
      _name("TrackerStatus"),_deadstraws(deadstraws) {}

    virtual ~TrackerStatus() {}

    std::string const& name() const { return _name; }
    void print( std::ostream& ) const;

    // find the status object(s) associated with a particular StrawId
    findMatches(StrawId const& sid, sstat_t& sstat) const {
      sstat.clear();
      for(auto istat= status_.begin(); istat != status_.end(); ++istat){
	if(istat->mask_.equal(sid,istat->sid_))sstat.append(istat->status_);
      }
    }

    // no signal expected from this straw
    bool noSignal(StrawId const& sid) {
      static StrawStatus::mask_type mask = StrawStatus::absent | StrawStatus::nowire | StrawStatus::noHV | StrawStatus::noLV | StrawStatus::nogas | StrawStatus::noPreamp;
      sstat_t matches;
      findMatches(sid,matches);
      bool retval(false);
      for(auto istat =matches.begin(); istat != matches.end(); istat++){
	retval |= istat->hasAnyProperty(mask);
      }
      return retval;
    }

    // whether to count in scattering calculation
    bool noMaterial(StrawStatus const& sstat) {
      static StrawStatus::mask_type mask = StrawStatus::absent;
      sstat_t matches;
      findMatches(sid,matches);
      bool retval(false);
      for(auto istat =matches.begin(); istat != matches.end(); istat++){
	retval |= istat->hasAnyProperty(mask);
      }
      return retval;
    }

  // suppress this straw
    bool suppress(StrawStatus const& sstat) {
      static StrawStatus::mask_type mask = StrawStatus::sparking | StrawStatus::noise | StrawStatus::pickup | StrawStatus::noread | StrawStatus::suppress;
      sstat_t matches;
      findMatches(sid,matches);
      bool retval(false);
      for(auto istat =matches.begin(); istat != matches.end(); istat++){
	retval |= istat->hasAnyProperty(mask);
      }
      return retval;
    }

  private:
    std::string _name;
    estat_t status_; // sparse list of tracker element status
  };

} // namespace mu2e

#endif /* TrackerConditions_TrackerStatus_hh */
