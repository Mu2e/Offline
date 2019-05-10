#ifndef TrackerConditions_FullReadoutStraw_hh
#define TrackerConditions_FullReadoutStraw_hh
//
// A list of straws that have full readout - no time window
// Initialized with FullReadoutStrawMaker
//

// Mu2e includes
#include "DataProducts/inc/StrawId.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "Mu2eInterfaces/inc/ProditionsEntity.hh"

// C++ includes
#include <set>
#include <memory>

namespace mu2e {

  class FullReadoutStraw : public ProditionsEntity {

  public:

    typedef std::shared_ptr<FullReadoutStraw> ptr_t;
    typedef std::shared_ptr<const FullReadoutStraw> cptr_t;
    typedef std::vector<StrawId> vec_t;

    FullReadoutStraw(): _name("FullReadoutStraw") {}
    FullReadoutStraw(vec_t const& straws):
      _name("FullReadoutStraw"),_straws(straws) {}

    virtual ~FullReadoutStraw() {}

    // Accessors; test if a straw is dead.  The position argument refers to 
    // distance from the straw center in mm
    bool isFull( StrawId id) const { return 
	std::find (_straws.begin(), _straws.end(), id)!=_straws.end();
    }

    std::string const& name() const { return _name; }
    void print( std::ostream& ) const;

  private:
    std::string _name;
    std::vector<StrawId> _straws; // sparse list of straws

  };

} // namespace mu2e

#endif /* TrackerConditions_FullReadoutStraw_hh */
