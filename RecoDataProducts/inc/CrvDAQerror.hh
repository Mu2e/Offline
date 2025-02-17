#ifndef RecoDataProducts_CrvDAQerror_hh
#define RecoDataProducts_CrvDAQerror_hh
//
//
// Contact person Ralf Ehrlich
//

#include "Offline/GeneralUtilities/inc/EnumToStringSparse.hh"
#include <vector>

namespace mu2e
{
  class CrvDAQerrorCodeDetail
  {
    public:

    enum enum_type{unknown=0, unableToGetDataBlock=1, invalidPacket=2, wrongSubsystemID=3, errorUnpackingStatusPacket=4, errorUnpackingCrvHits=5};
    static std::string const& typeName();
    static std::map<enum_type,std::string> const& names();
  };
  typedef EnumToStringSparse<CrvDAQerrorCodeDetail> CrvDAQerrorCode;

  class CrvDAQerror
  {
    public:

    CrvDAQerror() :
               _errorCode(), _subEvent(0), _dataBlock(0), _packetCount(0) {}

    CrvDAQerror(CrvDAQerrorCode::type errorCode, int subEvent, int dataBlock, int packetCount) :
               _errorCode(errorCode), _subEvent(subEvent), _dataBlock(dataBlock), _packetCount(packetCount) {}

    CrvDAQerrorCode::type GetErrorCode() const     {return _errorCode;}
    int                   GetSubEvent() const      {return _subEvent;}
    int                   GetDataBlock() const     {return _dataBlock;}
    int                   GetPacketCount() const   {return _packetCount;}

    private:

    CrvDAQerrorCode::type _errorCode;
    int                   _subEvent;
    int                   _dataBlock;
    int                   _packetCount;
  };
  typedef std::vector<mu2e::CrvDAQerror> CrvDAQerrorCollection;
}

#endif /* RecoDataProducts_CrvDAQerror_hh */
