#include "TrackerConditions/inc/TrackerDAQConditions.hh"

using namespace std;
namespace mu2e {

  StrawId TrackerDAQConditions::packetIdToStrawId(uint16_t packetId) const {
    //FIXME this is made up for now
    uint16_t dracChannel = packetId & 0xFF; 
    uint16_t rocNumber = (packetId & 0xFF00) >> 8;

    uint16_t straw = _DRACStrawMap[dracChannel];
    uint16_t upanel = _ROCPanelMap[rocNumber];
   
    uint16_t plane = (upanel / StrawId::_npanels);
    uint16_t panel = upanel - plane*StrawId::_npanels;

    return StrawId(plane,panel,straw);
  }

  void TrackerDAQConditions::print(std::ostream& os) const {
    os << endl << "TrackerDAQConditions parameters: "  << std::endl;
    os << "DRACStrawMap:" << endl;
    for (size_t i=0;i<_DRACStrawMap.size(); i++)
      os << "  " << _DRACStrawMap[i] << ", ";
    os << endl;
    os << "ROCPanelMap:" << endl;
    for (size_t i=0;i<_ROCPanelMap.size(); i++)
      os << "  " << _ROCPanelMap[i] << ", ";
    os << endl;
  }
}
