#ifndef RecoDataProducts_CrvStatus_hh
#define RecoDataProducts_CrvStatus_hh
//
//
// Contact person Ralf Ehrlich
//

#include "artdaq-core-mu2e/Overlays/Decoders/CRVDataDecoder.hh"

#include <vector>

namespace mu2e
{
  class CrvStatus
  {
    public:

      CrvStatus() {}

      //can't do a simple copy, because assigment and copy operators are implicitly deleted
      //for the DTC_DataHeaderPacket due the presence of a move constructor.
      CrvStatus(const DTCLib::DTC_DataHeaderPacket &dtcHeader, const DTCLib::DTC_SubEventHeader &subeventHeader)
        : _valid(),
        _linkID(static_cast<uint8_t>(dtcHeader.GetLinkID())),
        _eventWindowTag(dtcHeader.GetEventWindowTag().GetEventWindowTag(true)),
        _status(dtcHeader.GetStatus()),
        _dtcID(dtcHeader.GetID())
          //_subeventEventWindowTag(subeventHeader.GetEventWindowTag()),
          //_linkStatus(subeventHeader.GetLinkStatus()),
          //_linkLatency(subeventHeader.GetLinkLatency())
    {
      _subeventEventWindowTag = subeventHeader.event_tag_low +
        (static_cast<uint64_t>(subeventHeader.event_tag_high) << 32);

      if(dtcHeader.GetLinkID() == 0) {
        _linkStatus = subeventHeader.link0_status;
        _linkLatency = subeventHeader.link0_drp_rx_latency;
      } else if(dtcHeader.GetLinkID() == 1) {
        _linkStatus = subeventHeader.link1_status;
        _linkLatency = subeventHeader.link1_drp_rx_latency;
      } else if(dtcHeader.GetLinkID() == 2) {
        _linkStatus = subeventHeader.link2_status;
        _linkLatency = subeventHeader.link2_drp_rx_latency;
      } else if(dtcHeader.GetLinkID() == 3) {
        _linkStatus = subeventHeader.link3_status;
        _linkLatency = subeventHeader.link3_drp_rx_latency;
      } else if(dtcHeader.GetLinkID() == 4) {
        _linkStatus = subeventHeader.link4_status;
        _linkLatency = subeventHeader.link4_drp_rx_latency;
      } else if(dtcHeader.GetLinkID() == 5) {
        _linkStatus = subeventHeader.link5_status;
        _linkLatency = subeventHeader.link5_drp_rx_latency;
      }
      else {
        throw cet::exception("CRVStatus") << "Invalid link ID: " << dtcHeader.GetLinkID() << "\n";
      }
    }

      bool IsValid() const                         {return _valid;}
      uint8_t GetLinkID() const                    {return _linkID;}
      uint64_t GetEventWindowTag() const           {return _eventWindowTag;}
      uint8_t GetStatus() const                    {return _status;}
      uint8_t GetDTCID() const                     {return _dtcID;}
      uint64_t GetSubeventEventWindowTag() const   {return _subeventEventWindowTag;}
      uint8_t GetLinkStatus() const                {return _linkStatus;}
      uint16_t GetLinkLatency() const               {return _linkLatency;}
      std::vector<CRVDataDecoder::CRVROCStatusPacketFEBII> &GetROCHeader() {return _rocHeader;}
      const std::vector<CRVDataDecoder::CRVROCStatusPacketFEBII> &GetROCHeader() const {return _rocHeader;}
      bool HasROCHeader() const                    {return !_rocHeader.empty();}

      // --- MicroBunchStatus accessors (32-bit CRV ROC status word) ---
      // Bit layout from ROC firmware: https://github.com/Mu2e/CRVFirmware/tree/main/ROC/FPGA1
      // The ROC header is pushed after construction, so these read lazily
      // from _rocHeader[0].  Returns 0 / false when no ROC header is present.

      // Full 32-bit word
      uint32_t GetMicroBunchStatus() const {
        return _rocHeader.empty() ? 0 : _rocHeader[0].GetMicroBunchStatus();
      }

      // Bits [0:23] — per-port problem flags (3 groups of 8 ports)
      uint32_t GetPortFlags() const        {return GetMicroBunchStatus() & 0x00FFFFFFu;}

      // Bit 24 — any FEB micro-bunch number mismatch
      bool GetFEBuBMismatch() const        {return (GetMicroBunchStatus() >> 24) & 1;}
      // Bit 25 — any FEB buffer issue
      bool GetFEBBufferIssue() const       {return (GetMicroBunchStatus() >> 25) & 1;}
      // Bit 26 — any FEB overflow
      bool GetFEBOverflow() const          {return (GetMicroBunchStatus() >> 26) & 1;}
      // Bit 27 — group 1 issue (ports 0-7)
      bool GetGroup1Issue() const          {return (GetMicroBunchStatus() >> 27) & 1;}
      // Bit 28 — group 2 issue (ports 8-15)
      bool GetGroup2Issue() const          {return (GetMicroBunchStatus() >> 28) & 1;}
      // Bit 29 — group 3 issue (ports 16-23)
      bool GetGroup3Issue() const          {return (GetMicroBunchStatus() >> 29) & 1;}
      // Bit 30 — any group micro-bunch match error (includes ROC cross-link mismatch)
      bool GetuBMatchError() const         {return (GetMicroBunchStatus() >> 30) & 1;}
      // Bit 31 — any group truncation
      bool GetTruncation() const           {return (GetMicroBunchStatus() >> 31) & 1;}

    private:

      // Data Header Packet
      bool                        _valid{false};
      uint8_t                     _linkID;
      uint64_t                    _eventWindowTag;
      uint8_t                     _status{0};
      uint8_t                     _dtcID{0};
      // Subevent Packet
      uint64_t                     _subeventEventWindowTag;
      uint8_t                     _linkStatus;
      uint16_t                    _linkLatency;

      std::vector<CRVDataDecoder::CRVROCStatusPacketFEBII>  _rocHeader;   //not every DTC link has a ROC attached.
      //actually wanted to use std::optional,
      //but this causes problems with the dictionary.
  };
  typedef std::vector<CrvStatus> CrvStatusCollection;  //one entry per block
}

#endif /* RecoDataProducts_CrvStatus_hh */
