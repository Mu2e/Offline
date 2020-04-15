#ifndef TrackerConditions_TrackerDAQConditions_hh
#define TrackerConditions_TrackerDAQConditions_hh

//
// TrackerDAQConditions stores channel maps for packet parsing
//

// C++ includes
#include <vector>

// Mu2e includes
#include "DataProducts/inc/StrawId.hh"
#include "DataProducts/inc/TrackerPacketId.hh"
#include "Mu2eInterfaces/inc/ProditionsEntity.hh"
#include "fhiclcpp/ParameterSet.h"

namespace mu2e {

  class TrackerDAQConditions : virtual public ProditionsEntity {
  public:

    typedef std::shared_ptr<TrackerDAQConditions> ptr_t;
    typedef std::shared_ptr<const TrackerDAQConditions> cptr_t;

    TrackerDAQConditions():_name("TrackerDAQConditions") {}

    // construct with constants, then some values are computed and filled below
    TrackerDAQConditions(std::vector<uint16_t> DRACStrawMap, std::vector<uint16_t> ROCPanelMap) :
      _name("TrackerDAQConditions"),
      _DRACStrawMap(DRACStrawMap), _ROCPanelMap(ROCPanelMap){}

    virtual ~TrackerDAQConditions() {}
    
    std::string const& name() const { return _name; }
    void print(std::ostream& os) const;

    StrawId packetIdToStrawId(TrackerPacketId const& packetId) const;
    TrackerPacketId strawIdToPacketId(StrawId const& strawId) const;

    // all of these must be called to fill this object
    void setDRACStrawMap(std::vector<uint16_t> DRACStrawMap) { _DRACStrawMap = DRACStrawMap; }
    void setROCPanelMap(std::vector<uint16_t> ROCPanelMap) { _ROCPanelMap = ROCPanelMap; }

  private:
    
    std::string _name;

    std::vector<uint16_t> _DRACStrawMap;
    std::vector<uint16_t> _ROCPanelMap;

  };
  
}

#endif

