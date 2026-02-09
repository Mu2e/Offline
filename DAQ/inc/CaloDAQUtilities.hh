#ifndef DAQ_CaloDAQUtilities_hh
#define DAQ_CaloDAQUtilities_hh

//
//
//

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "artdaq-core-mu2e/Overlays/Decoders/CalorimeterDataDecoder.hh"
#include "artdaq-core-mu2e/Overlays/FragmentType.hh"
#include <artdaq-core/Data/ContainerFragment.hh>
#include <artdaq-core/Data/Fragment.hh>

#include "Offline/DataProducts/inc/CaloConst.hh"

#include <string>
#include <vector>

namespace mu2e {

class CaloDAQUtilities {

public:
  CaloDAQUtilities(std::string ModuleName);

  enum CaloHitError {
    Good = 0,
    BeginMarker = 1,
    LastSampleMarker = 2,
    WaveformSize = 3,
    NumberOfSamples = 4,
    MaxSampleIndex = 5,
    BoardID = 6,
    ChannelID = 7,
    ErrorFlags = 8
  };

  std::string getCaloHitErrorName(CaloHitError error) {
    switch (error) {
    case Good:
      return "Good";
    case BeginMarker:
      return "BeginMarker";
    case LastSampleMarker:
      return "LastSampleMarker";
    case WaveformSize:
      return "WaveformSize";
    case NumberOfSamples:
      return "NumberOfSamples";
    case MaxSampleIndex:
      return "MaxSampleIndex";
    case BoardID:
      return "BoardID";
    case ChannelID:
      return "ChannelID";
    case ErrorFlags:
      return "ErrorFlags";
    default:
      return "Unknown";
    }
  }


  CaloHitError isHitGood(std::pair<CalorimeterDataDecoder::CalorimeterHitDataPacket, std::vector<uint16_t>> const& Hit) {
    if (Hit.first.Reserved1 != 0xAAA)
      return CaloHitError::BeginMarker;
    if (Hit.first.ErrorFlags != 0)
      return CaloHitError::ErrorFlags;
    if (Hit.second.size() == 0)
      return CaloHitError::WaveformSize;
    if (Hit.first.NumberOfSamples != Hit.second.size())
      return CaloHitError::NumberOfSamples;
    if (Hit.first.IndexOfMaxDigitizerSample >= Hit.second.size())
      return CaloHitError::MaxSampleIndex;
    if (Hit.first.BoardID >= CaloConst::_nDIRAC)
      return CaloHitError::BoardID;
    if (Hit.first.ChannelID >= CaloConst::_nChPerDIRAC)
      return CaloHitError::ChannelID;
    return CaloHitError::Good;
  }

  CaloHitError isHitGood(std::pair<CalorimeterDataDecoder::CalorimeterHitDataPacket, uint16_t> const& Hit) {
    if (Hit.first.Reserved1 != 0xAAA)
      return CaloHitError::BeginMarker;
    if (Hit.first.ErrorFlags != 0)
      return CaloHitError::ErrorFlags;
    if (Hit.second == 0)
      return CaloHitError::WaveformSize;
    if (Hit.first.BoardID >= CaloConst::_nDIRAC)
      return CaloHitError::BoardID;
    if (Hit.first.ChannelID >= CaloConst::_nChPerDIRAC)
      return CaloHitError::ChannelID;
    return CaloHitError::Good;
  }

  CaloHitError isHitGood(std::pair<CalorimeterDataDecoder::CalorimeterHitTestDataPacket, std::vector<uint16_t>> const& Hit) {
    if (Hit.first.BeginMarker != 0xAAA)
      return CaloHitError::BeginMarker;
    if (Hit.first.LastSampleMarker == 0)
      return CaloHitError::LastSampleMarker;
    if (Hit.second.size() == 0)
      return CaloHitError::WaveformSize;
    if (Hit.first.NumberOfSamples != Hit.second.size())
      return CaloHitError::NumberOfSamples;
    if (Hit.first.IndexOfMaxDigitizerSample >= Hit.second.size())
      return CaloHitError::MaxSampleIndex;
    if (Hit.first.BoardID >= CaloConst::_nDIRAC)
      return CaloHitError::BoardID;
    if (Hit.first.ChannelID >= CaloConst::_nChPerDIRAC)
      return CaloHitError::ChannelID;
    return CaloHitError::Good;
  }

  CaloHitError isHitGood(std::pair<CalorimeterDataDecoder::CalorimeterHitTestDataPacket, uint16_t> const& Hit) {
    if (Hit.first.BeginMarker != 0xAAA)
      return CaloHitError::BeginMarker;
    if (Hit.first.LastSampleMarker == 0)
      return CaloHitError::LastSampleMarker;
    if (Hit.second == 0)
      return CaloHitError::WaveformSize;
    if (Hit.first.BoardID >= CaloConst::_nDIRAC)
      return CaloHitError::BoardID;
    if (Hit.first.ChannelID >= CaloConst::_nChPerDIRAC)
      return CaloHitError::ChannelID;
    return CaloHitError::Good;
  }

  void printCaloFragmentInfo(CalorimeterDataDecoder const& Frag);

  void printCaloFragmentHeader(std::shared_ptr<DTCLib::DTC_DataHeaderPacket> Header);

  void printCaloPulse(CalorimeterDataDecoder::CalorimeterHitDataPacket const& Hit);
  void printCaloPulse(CalorimeterDataDecoder::CalorimeterHitTestDataPacket const& Hit);

  void printWaveform(std::vector<uint16_t> const& Pulse);


  // Function to get art fragments from event
  artdaq::Fragments getFragments(art::Event& event) {

    artdaq::Fragments fragments;
    artdaq::FragmentPtrs containerFragments;

    std::vector<art::Handle<artdaq::Fragments>> fragmentHandles;
    fragmentHandles = event.getMany<std::vector<artdaq::Fragment>>();
    for (const auto& handle : fragmentHandles) {
      if (!handle.isValid() || handle->empty()) {
        continue;
      }

      if (handle->front().type() == artdaq::Fragment::ContainerFragmentType) {
        for (const auto& cont : *handle) {
          artdaq::ContainerFragment contf(cont);
          if (contf.fragment_type() != mu2e::FragmentType::DTCEVT) {
            break;
          }

          for (size_t ii = 0; ii < contf.block_count(); ++ii) {
            containerFragments.push_back(contf[ii]);
            fragments.push_back(*containerFragments.back());
          }
        }
      } else {
        if (handle->front().type() == mu2e::FragmentType::DTCEVT) {
          for (auto frag : *handle) {
            fragments.emplace_back(frag);
          }
        }
      }
    }
    return fragments;
  }

private:
  std::string moduleName_;
};
} // namespace mu2e

#endif /* DAQ_CaloDAQUtilities_hh */
