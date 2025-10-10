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
    CrvStatus(const DTCLib::DTC_DataHeaderPacket &dtcHeader)
	: _valid(),
	  _linkID(static_cast<uint8_t>(dtcHeader.GetLinkID())), 
          _eventWindowTag(dtcHeader.GetEventWindowTag().GetEventWindowTag(true)),
          _status(dtcHeader.GetStatus()),
          _dtcID(dtcHeader.GetID()) 
	  {}

    bool IsValid() const                         {return _valid;}
    uint8_t GetLinkID() const                    {return _linkID;}
    uint64_t GetEventWindowTag() const           {return _eventWindowTag;}
    uint8_t GetStatus() const                    {return _status;}
    uint8_t GetDTCID() const                     {return _dtcID;}
    std::vector<CRVDataDecoder::CRVROCStatusPacket> &GetROCHeader() {return _rocHeader;}

    private:

    bool                        _valid{false};
    uint8_t                     _linkID;
    uint64_t                    _eventWindowTag;
    uint8_t                     _status{0};
    uint8_t                     _dtcID{0};
    std::vector<CRVDataDecoder::CRVROCStatusPacket>  _rocHeader;   //not every DTC link has a ROC attached.
                                                                   //actually wanted to use std::optional, 
								   //but this causes problems with the dictionary.
  };
  typedef std::vector<CrvStatus> CrvStatusCollection;  //one entry per block
}

#endif /* RecoDataProducts_CrvStatus_hh */

