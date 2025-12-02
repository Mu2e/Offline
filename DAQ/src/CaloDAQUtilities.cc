
#include "artdaq-core-mu2e/Overlays/Decoders/CalorimeterDataDecoder.hh"
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

// clang-format off
void CaloDAQUtilities::printCaloFragmentHeader(std::shared_ptr<DTCLib::DTC_DataHeaderPacket> Header) {
  std::cout << "timestamp: "           << static_cast<int>(Header->GetEventWindowTag().GetEventWindowTag(true)) << std::endl;
  std::cout << "Header->SubsystemID: " << static_cast<int>(Header->GetSubsystemID()) << std::endl;
  std::cout << "dtcID: "               << static_cast<int>(Header->GetID()) << std::endl;
  std::cout << "rocID: "               << static_cast<int>(Header->GetLinkID()) << std::endl;
  std::cout << "packetCount: "         << static_cast<int>(Header->GetPacketCount()) << std::endl;
  std::cout << "EVB mode: "            << static_cast<int>(Header->GetEVBMode()) << std::endl;
}

void CaloDAQUtilities::printCaloPulse(CalorimeterDataDecoder::CalorimeterHitDataPacket const& Hit) {
  std::cout << "[CaloDAQUtilities] \tBoardID                    " << (int)Hit.BoardID << std::endl;
  std::cout << "[CaloDAQUtilities] \tDetectorID                 " << (int)Hit.DetectorID<< std::endl;
  std::cout << "[CaloDAQUtilities] \tChannelID                  " << (int)Hit.ChannelID << std::endl;
  std::cout << "[CaloDAQUtilities] \tTime                       " << (int)Hit.Time << std::endl;
  std::cout << "[CaloDAQUtilities] \tInPayloadEventWindowTag    " << (int)Hit.InPayloadEventWindowTag << std::endl;
  std::cout << "[CaloDAQUtilities] \tBaseline                   " << (int)Hit.Baseline << std::endl;
  std::cout << "[CaloDAQUtilities] \tIndexOfMaxDigitizerSample  " << (int)Hit.IndexOfMaxDigitizerSample << std::endl;
  std::cout << "[CaloDAQUtilities] \tErrorFlags                 " << (int)Hit.ErrorFlags<< std::endl;
  std::cout << "[CaloDAQUtilities] \tNumberOfSamples            " << (int)Hit.NumberOfSamples<< std::endl;
}

void CaloDAQUtilities::printCaloPulse(CalorimeterDataDecoder::CalorimeterHitTestDataPacket const& Hit) {
  std::cout << "[CaloDAQUtilities] \tBoardID      " << (int)Hit.BoardID << std::endl;
  std::cout << "[CaloDAQUtilities] \tChNumber     " << (int)Hit.ChannelID << std::endl;
  std::cout << "[CaloDAQUtilities] \tEWT          " << (int)Hit.InPayloadEventWindowTag<< std::endl;
  std::cout << "[CaloDAQUtilities] \tErrorFlags   " << (int)Hit.ErrorFlags << std::endl;
  std::cout << "[CaloDAQUtilities] \tTime         " << (int)Hit.Time << std::endl;
  std::cout << "[CaloDAQUtilities] \tNSamples     " << (int)Hit.NumberOfSamples << std::endl;
  std::cout << "[CaloDAQUtilities] \tIndexMax     " << (int)Hit.IndexOfMaxDigitizerSample<< std::endl;
}
// clang-format on

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

} // namespace mu2e
