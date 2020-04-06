

#include "TrackerConditions/inc/TrackerDAQConditionsMaker.hh"
#include "DataProducts/inc/StrawId.hh"

#include <vector>

using namespace std;

namespace mu2e {

  TrackerDAQConditions::ptr_t TrackerDAQConditionsMaker::fromFcl() {

    // creat this at the beginning since it must be used,
    // partially constructed, to complete the construction
    auto ptr = std::make_shared<TrackerDAQConditions>(_config.DRACStrawMap(), _config.ROCPanelMap());

    return ptr;

  } // end fromFcl

  TrackerDAQConditions::ptr_t TrackerDAQConditionsMaker::fromDb(
      TrkDRACtoStraw::cptr_t tdts,
      TrkROCtoPanel::cptr_t trtp ) {

    // initially fill from fcl to get all the constants
    auto ptr = fromFcl();

    // now overwrite with db values
    vector<uint16_t> dracstrawmap(StrawId::_nstraws);
    for (size_t i=0;i<StrawId::_nstraws;i++){
      dracstrawmap[i] = tdts->rowAt(i).straw();
    }
    ptr->setDRACStrawMap(dracstrawmap);

    vector<uint16_t> rocpanelmap(StrawId::_nupanels);
    for (size_t i=0;i<StrawId::_nupanels;i++){
      rocpanelmap[i] = trtp->rowAt(i).panel();
    }
    ptr->setROCPanelMap(rocpanelmap);
    
    return ptr;

  } // end fromDb

}
