
#include "artdaq-core-mu2e/Data/CalorimeterDataDecoder.hh"
#include "artdaq-core-mu2e/Overlays/FragmentType.hh"
#include <artdaq-core/Data/Fragment.hh>

#include "Offline/DAQ/inc/CaloDAQUtilities.hh"
#include <iostream>

namespace mu2e {
CaloDAQUtilities::CaloDAQUtilities(std::string ModuleName) : moduleName_(ModuleName) {}

void CaloDAQUtilities::printCaloFragmentInfo(CalorimeterDataDecoder const& cc) {
  std::cout << std::endl;
  std::cout << "ArtFragmentReader: ";
  std::cout << "\tBlock Count: " << std::dec << cc.block_count() << std::endl;
  std::cout << std::endl;
  std::cout << "\t"
            << "====== Example Block Sizes ======" << std::endl;
  for (size_t i = 0; i < 10; i++) {
    if (i < cc.block_count()) {
      std::cout << "\t" << i << "\t" << cc.blockSizeBytes(i) << std::endl;
    }
  }
  std::cout << "\t"
            << "=========================" << std::endl;
}

void CaloDAQUtilities::printCaloFragmentHeader(
    std::shared_ptr<DTCLib::DTC_DataHeaderPacket> Header) {

  std::cout << "timestamp: "
            << static_cast<int>(Header->GetEventWindowTag().GetEventWindowTag(true)) << std::endl;
  std::cout << "Header->SubsystemID: " << static_cast<int>(Header->GetSubsystemID()) << std::endl;
  std::cout << "dtcID: " << static_cast<int>(Header->GetID()) << std::endl;
  std::cout << "rocID: " << static_cast<int>(Header->GetLinkID()) << std::endl;
  std::cout << "packetCount: " << static_cast<int>(Header->GetPacketCount()) << std::endl;
  std::cout << "EVB mode: " << static_cast<int>(Header->GetEVBMode()) << std::endl;

  std::cout << std::endl;
}

void CaloDAQUtilities::printCaloPulse(CalorimeterDataDecoder::CalorimeterHitDataPacket const& Hit) {
  std::cout << "[CaloDAQUtilities] \tChNumber   " << (int)Hit.ChannelNumber << std::endl;
  std::cout << "[CaloDAQUtilities] \tDIRACA     " << (int)Hit.DIRACA << std::endl;
  std::cout << "[CaloDAQUtilities] \tDIRACB     " << (int)Hit.DIRACB << std::endl;
  std::cout << "[CaloDAQUtilities] \tErrorFlags " << (int)Hit.ErrorFlags << std::endl;
  std::cout << "[CaloDAQUtilities] \tTime              " << (int)Hit.Time << std::endl;
  std::cout << "[CaloDAQUtilities] \tNSamples   " << (int)Hit.NumberOfSamples << std::endl;
  std::cout << "[CaloDAQUtilities] \tIndexMax   " << (int)Hit.IndexOfMaxDigitizerSample
            << std::endl;
}

void CaloDAQUtilities::printCaloPulse(
    CalorimeterDataDecoder::CalorimeterHitTestDataPacket const& Hit) {
  std::cout << "[CaloDAQUtilities] \tBoardID      " << (int)Hit.BoardID << std::endl;
  std::cout << "[CaloDAQUtilities] \tChNumber     " << (int)Hit.ChannelID << std::endl;
  std::cout << "[CaloDAQUtilities] \tEWT          " << (int)Hit.InPayloadEventWindowTag
            << std::endl;
  std::cout << "[CaloDAQUtilities] \tErrorFlags   " << (int)Hit.ErrorFlags << std::endl;
  std::cout << "[CaloDAQUtilities] \tTime         " << (int)Hit.Time << std::endl;
  std::cout << "[CaloDAQUtilities] \tNSamples     " << (int)Hit.NumberOfSamples << std::endl;
  std::cout << "[CaloDAQUtilities] \tIndexMax     " << (int)Hit.IndexOfMaxDigitizerSample
            << std::endl;
}

void CaloDAQUtilities::printWaveform(std::vector<uint16_t> const& Pulse) {
  std::cout << "Waveform: {";
  for (size_t i = 0; i < Pulse.size(); i++) {
    std::cout << Pulse[i];
    if (i < Pulse.size() - 1) {
      std::cout << ",";
    }
  }
  std::cout << "}" << std::endl;
}

void CaloDAQUtilities::printAllHitInfo(int CrystalID, int SiPMID,
                                       std::shared_ptr<DTCLib::DTC_DataHeaderPacket> Header,
                                       CalorimeterDataDecoder::CalorimeterHitDataPacket const& Hit,
                                       uint16_t PulseMax) {

  std::cout << "Crystal ID: " << CrystalID << std::endl;
  std::cout << "SiPM ID: " << SiPMID << std::endl;
  std::cout << "Time: " << (int)Hit.Time << std::endl;
  std::cout << "NumSamples: " << (int)Hit.NumberOfSamples << std::endl;
  // printWaveform(Pulse);

  // Text format: timestamp crystalID roID time nsamples samples...
  // Example: 1 201 402 660 18 0 0 0 0 1 17 51 81 91 83 68 60 58 52 42 33 23 16
  std::cout << "GREPMECAL: " << Header->GetEventWindowTag().GetEventWindowTag(true) << " ";
  std::cout << CrystalID << " ";
  std::cout << SiPMID << " ";
  std::cout << Hit.Time << " ";
  std::cout << PulseMax << " ";
  // std::cout << Pulse.size()   << " ";
  // for (size_t i = 0; i < Pulse.size(); i++) {
  //   std::cout << Pulse[i];
  //   if (i < Pulse.size() - 1) {
  //         std::cout << " ";
  //   }
  // }
  std::cout << std::endl;
}

} // namespace mu2e
