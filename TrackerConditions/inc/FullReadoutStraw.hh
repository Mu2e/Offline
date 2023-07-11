#ifndef TrackerConditions_FullReadoutStraw_hh
#define TrackerConditions_FullReadoutStraw_hh
//
// A list of straws that have full readout - no time window
// Initialized with FullReadoutStrawMaker
//

// Mu2e includes
#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/TrackerGeom/inc/Straw.hh"
#include "Offline/Mu2eInterfaces/inc/ProditionsEntity.hh"

// C++ includes
#include <set>
#include <memory>

namespace mu2e {

  class FullReadoutStraw : public ProditionsEntity {

    public:

      typedef std::shared_ptr<FullReadoutStraw> ptr_t;
      typedef std::shared_ptr<const FullReadoutStraw> cptr_t;
      typedef std::vector<StrawId> vec_t;
      constexpr static const char* cxname = {"FullReadoutStraw"};

      FullReadoutStraw(vec_t const& straws):
        ProditionsEntity(cxname),_straws(straws) {}

      virtual ~FullReadoutStraw() {}

      // Accessors; test if a straw is dead.  The position argument refers to
      // distance from the straw center in mm
      bool isFull( StrawId id) const { return
        std::find (_straws.begin(), _straws.end(), id)!=_straws.end();
      }

      void print( std::ostream& ) const;

    private:
      std::vector<StrawId> _straws; // sparse list of straws

  };

} // namespace mu2e

#endif /* TrackerConditions_FullReadoutStraw_hh */
