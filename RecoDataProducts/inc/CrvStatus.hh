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
    uint8_t GetLinkLatency() const               {return _linkLatency;}
    std::vector<CRVDataDecoder::CRVROCStatusPacketFEBII> &GetROCHeader() {return _rocHeader;}

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
    uint8_t                     _linkLatency;

    std::vector<CRVDataDecoder::CRVROCStatusPacketFEBII>  _rocHeader;   //not every DTC link has a ROC attached.
                                                                        //actually wanted to use std::optional, 
							                //but this causes problems with the dictionary.
  };
  typedef std::vector<CrvStatus> CrvStatusCollection;  //one entry per block
}

#endif /* RecoDataProducts_CrvStatus_hh */

