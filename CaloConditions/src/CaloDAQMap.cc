#include "Offline/CaloConditions/inc/CaloDAQMap.hh"
#include <stdexcept>

using namespace std;
namespace mu2e {

  uint16_t CaloDAQMap::packetIdTocaloRoId(uint16_t packetId) const {
    //
    // FIXME this is made up for now FF-mask for channel Number and BoardNumber
    // DiracChannel (0-19) OK in FF , DiracNumber (0-136) ok in FF
    // No INC-files to descrive Board numbering yet. First Check V09 mapping
    //
    /* version 09
    uint16_t DiracChannel = packetId & 0xFF;
    uint16_t DiracNumber = (packetId & 0xFF00) >> 8; */

    uint16_t DiracNumber  = packetId & 0xFF;              // first 8 bits
    uint16_t DiracChannel = (packetId & 0x1F00) >> 8;     // second 5 bits
    uint16_t DetType      = (packetId & 0x7000) >> 13;    // last 3 bits: 0=Calo,1=Caphri,2=Pin

    if (DiracChannel >= 19 || DiracNumber >= 136){
      throw cet::exception("BADINPUTS")<<"CaloDAQMap::packetIdTocaloRoId : packetId DiracChannel/DiracNumber out of range" << std::endl;
    }

    if( DetType <=1 ){
      uint16_t DiracPoi = DiracNumber*20 + DiracChannel;
      uint16_t roId  = _DIRAC2CaloMap[DiracPoi];
      if( roId < 0 ){
        throw cet::exception("BadRoId")<<"CaloDAQMap::packetIdTocaloRoId : roID for Calo empty channels or Calo PINs" << std::endl;
      }
      //printf(" DIRAC # %d Chan # %d -> RoId %d \n",DiracChannel,DiracNumber,roId);
      return roId;
    }else{
      throw cet::exception("BadRoId")<<"CaloDAQMap::packetIdTocaloRoId : DetType not implemented" << std::endl;
    }
  }

  uint16_t CaloDAQMap::caloRoIdToPacketId(uint16_t caloRoId) const {

    if ( caloRoId >= 674*4  ){
      throw cet::exception("BADINPUTS")<<"CaloDAQMap::caloRoIdToPacketId : caloRoId out of range" << std::endl;
    }

    uint16_t DiracPoi  = _Calo2DIRACMap[caloRoId];
    if ( DiracPoi >= 136*20 ){
      throw cet::exception("BADINPUTS")<<"CaloDAQMap::caloRoIdToPacketId : return DiracPoi  out of range" << std::endl;
    }

    uint16_t DiracNum = DiracPoi/20;
    uint16_t DiracChannel = DiracPoi -20*DiracNum;

    // printf(" RoID %d --> DIRAC Poi %d DIRAC # %d Chan # %d \n",
    //           caloRoId,DiracPoi,DiracNum,DiracChannel);

    //
    // Now from CaloRoId decide DetType to discriminate btw Calo/Caphri
    // Leave PIN-diodes for another round
    //
    uint16_t DetType =0;
    if ( caloRoId == 624 || caloRoId == 625) DetType=1; // FIX ME As soon as I know it when and how
    uint16_t PacketId;
    PacketId = (DiracNum & 0xFF) | ((DiracChannel & 0x1F) << 8)| ((DetType & 0x7) << 13);
    return PacketId;
  }

  void CaloDAQMap::print(std::ostream& os) const {
    os << endl << "CaloDAQMap parameters: "  << std::endl;
    os << "DIRAC2CaloMap:" << endl;
    for (size_t i=0;i<_DIRAC2CaloMap.size(); i++)
      os << "  " << _DIRAC2CaloMap[i] << ", ";
    os << endl;
    os << "Calo2DIRACMap:" << endl;
    for (size_t i=0;i<_Calo2DIRACMap.size(); i++)
      os << "  " << _Calo2DIRACMap[i] << ", ";
    os << endl;
  }
}
