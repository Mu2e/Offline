#ifndef DAQ_CaloDAQUtilities_hh
#define DAQ_CaloDAQUtilities_hh

//
//
//

#include "artdaq-core-mu2e/Data/CalorimeterDataDecoder.hh"
#include "artdaq-core-mu2e/Overlays/FragmentType.hh"
#include <artdaq-core/Data/Fragment.hh>

#include <string>

namespace mu2e {

  class CaloDAQUtilities{

  public:
    CaloDAQUtilities(std::string ModuleName);

    uint16_t  getCrystalID(CalorimeterDataDecoder::CalorimeterHitDataPacket const& Hit){ return Hit.DIRACB & 0x0FFF;}
    uint16_t  getSiPMID   (CalorimeterDataDecoder::CalorimeterHitDataPacket const& Hit){
      uint16_t  crystalID  = getCrystalID(Hit);
      uint16_t  sipmID     = Hit.DIRACB >> 12;
      return (crystalID * 2 + sipmID);
    }

    void   printCaloFragmentInfo(CalorimeterDataDecoder const& Frag);

    void   printCaloFragmentHeader(std::shared_ptr<DTCLib::DTC_DataHeaderPacket> Header);

    void   printCaloPulse(CalorimeterDataDecoder::CalorimeterHitDataPacket const& Hit);

    void   printWaveform(std::vector<uint16_t> const& Pulse);

    void   printAllHitInfo(int CrystalID, int SiPMID, std::shared_ptr<DTCLib::DTC_DataHeaderPacket> Header, CalorimeterDataDecoder::CalorimeterHitReadoutPacket const& Hit, uint16_t PulseMax);

  private:
    std::string  moduleName_;
  };
}

#endif /* DAQ_CaloDAQUtilities_hh */
