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

    enum CaloHitError {
      Good = 0,
      BeginMarker = 1,
      LastSampleMarker = 2,
      WaveformSize = 3,
      NumberOfSamples = 4,
      MaxSampleIndex = 5,
      BoardID = 6,
      ChannelID = 7
    };

    std::string getCaloHitErrorName(CaloHitError error) {
      switch (error) {
        case Good: return "Good";
        case BeginMarker: return "BeginMarker";
        case LastSampleMarker: return "LastSampleMarker";
        case WaveformSize: return "WaveformSize";
        case NumberOfSamples: return "NumberOfSamples";
        case MaxSampleIndex: return "MaxSampleIndex";
        case BoardID: return "BoardID";
        case ChannelID: return "ChannelID";
        default: return "Unknown";
      }
    }

    CaloHitError isHitGood(std::pair<CalorimeterDataDecoder::CalorimeterHitTestDataPacket, std::vector<uint16_t>> const& Hit){
      if (Hit.first.BeginMarker != 0xAAA) return CaloHitError::BeginMarker;
      if (Hit.first.LastSampleMarker == 0) return CaloHitError::LastSampleMarker;
      if (Hit.second.size() == 0) return CaloHitError::WaveformSize;
      if (Hit.first.NumberOfSamples != Hit.second.size()) return CaloHitError::NumberOfSamples;
      if (Hit.first.IndexOfMaxDigitizerSample >= Hit.second.size()) return CaloHitError::NumberOfSamples;
      if (Hit.first.BoardID < 0 || Hit.first.BoardID >= 160) return CaloHitError::BoardID;
      if (Hit.first.ChannelID < 0 || Hit.first.ChannelID >= 20) return CaloHitError::ChannelID;
      return CaloHitError::Good;
    }

    void   printCaloFragmentInfo(CalorimeterDataDecoder const& Frag);

    void   printCaloFragmentHeader(std::shared_ptr<DTCLib::DTC_DataHeaderPacket> Header);

    void   printCaloPulse(CalorimeterDataDecoder::CalorimeterHitDataPacket const& Hit);
    void   printCaloPulse(CalorimeterDataDecoder::CalorimeterHitTestDataPacket const& Hit);

    void   printWaveform(std::vector<uint16_t> const& Pulse);

    void   printAllHitInfo(int CrystalID, int SiPMID, std::shared_ptr<DTCLib::DTC_DataHeaderPacket> Header, CalorimeterDataDecoder::CalorimeterHitDataPacket const& Hit, uint16_t PulseMax);

  private:
    std::string  moduleName_;
  };
}

#endif /* DAQ_CaloDAQUtilities_hh */
