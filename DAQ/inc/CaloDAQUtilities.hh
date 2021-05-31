#ifndef DAQ_CaloDAQUtilities_hh
#define DAQ_CaloDAQUtilities_hh

//
//
//

#include "mu2e-artdaq-core/Overlays/CalorimeterFragment.hh"
#include "mu2e-artdaq-core/Overlays/FragmentType.hh"
#include <artdaq-core/Data/Fragment.hh>

#include <string>

namespace mu2e {

  class CaloDAQUtilities{

  public:
    CaloDAQUtilities(std::string ModuleName);

    uint16_t  getCrystalID(CalorimeterFragment::CalorimeterHitReadoutPacket& Hit){ return Hit.DIRACB & 0x0FFF;}
    uint16_t  getSiPMID   (CalorimeterFragment::CalorimeterHitReadoutPacket& Hit){ 
      uint16_t  crystalID  = getCrystalID(Hit);
      uint16_t  sipmID     = Hit.DIRACB >> 12;
      return (crystalID * 2 + sipmID);
    }

    void   printCaloFragmentInfo(const artdaq::Fragment& f, CalorimeterFragment& Frag);

    void   printCaloFragmentHeader(DTCLib::DTC_DataHeaderPacket &Header);
  
    void   printCaloPulse(CalorimeterFragment::CalorimeterHitReadoutPacket& Hit);

    void   printWaveform(std::vector<uint16_t>& Pulse);

    void   printAllHitInfo(int CrystalID, int SiPMID, DTCLib::DTC_DataHeaderPacket &Header, CalorimeterFragment::CalorimeterHitReadoutPacket& Hit, uint16_t& PulseMax);
    
  private:
    std::string  moduleName_;
  };
}

#endif /* DAQ_CaloDAQUtilities_hh */
