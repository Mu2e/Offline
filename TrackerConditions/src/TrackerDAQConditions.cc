#include "TrackerConditions/inc/TrackerDAQConditions.hh"
#include <stdexcept>

using namespace std;
namespace mu2e {

  StrawId TrackerDAQConditions::packetIdToStrawId(TrackerPacketId const& packetId) const {
    uint16_t dracChannel = packetId.channel(); 
    uint16_t rocNumber = packetId.roc();

    if (dracChannel >= _DRACStrawMap.size() || rocNumber >= _ROCPanelMap.size()){
      throw cet::exception("BADINPUTS")<<"TrackerDAQConditions::packetIdToStrawId : packetId dracChannel/rocNumber out of range" << std::endl;
    }

    uint16_t straw = _DRACStrawMap[dracChannel];
    uint16_t upanel = _ROCPanelMap[rocNumber];
   
    uint16_t plane = (upanel / StrawId::_npanels);
    uint16_t panel = upanel - plane*StrawId::_npanels;

    return StrawId(plane,panel,straw);
  }
  
  TrackerPacketId TrackerDAQConditions::strawIdToPacketId(StrawId const& strawId) const {
    auto iterStraw = std::find(_DRACStrawMap.begin(),_DRACStrawMap.end(),strawId.getStraw());
    if (iterStraw == _DRACStrawMap.end())
      throw cet::exception("BADINPUTS")<<"TrackerDAQConditions::strawIdToPacketId : strawId out of range" << std::endl;
    auto iterPanel = std::find(_ROCPanelMap.begin(),_ROCPanelMap.end(),strawId.uniquePanel());
    if (iterPanel == _ROCPanelMap.end())
      throw cet::exception("BADINPUTS")<<"TrackerDAQConditions::strawIdToPacketId : strawId out of range" << std::endl;

    return TrackerPacketId(*iterStraw, *iterPanel);
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
