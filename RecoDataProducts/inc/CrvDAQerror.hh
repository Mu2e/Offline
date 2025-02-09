#ifndef RecoDataProducts_CrvDAQerror_hh
#define RecoDataProducts_CrvDAQerror_hh
//
//
// Contact person Ralf Ehrlich
//

#include <vector>

namespace mu2e
{
  class CrvDAQerror
  {
    public:

    enum errorCode{unableToGetDataBlock=0, invalidPacket=1, wrongSubsystemID=2, errorUnpackingStatusPacket=3, errorUnpackingCrvHits=4};

    CrvDAQerror() {}

    CrvDAQerror(int errorCode, int subEvent, int dataBlock, int packetCount) :
               _errorCode(errorCode), _subEvent(subEvent), _dataBlock(dataBlock), _packetCount(packetCount) {}

    int  GetErrorCode() const     {return _errorCode;}
    int  GetSubEvent() const      {return _subEvent;}
    int  GetDataBlock() const     {return _dataBlock;}
    int  GetPacketCount() const   {return _packetCount;}

    private:

    int    _errorCode{0};
    int    _subEvent{0};
    int    _dataBlock{0};
    int    _packetCount{0};
  };
  typedef std::vector<mu2e::CrvDAQerror> CrvDAQerrorCollection;
}

#endif /* RecoDataProducts_CrvDAQerror_hh */
