#include <type_traits>
#include <utility>

#include "Offline/RecoDataProducts/inc/CrvDAQerror.hh"

namespace mu2e
{
  std::string const& CrvDAQerrorCodeDetail::typeName()
  {
    static const std::string type("CrvDAQerrorCode");
    return type;
  }

  static const std::map<CrvDAQerrorCodeDetail::enum_type,std::string> nam
  {
    {CrvDAQerrorCodeDetail::unknown                    , "unknown" },
    {CrvDAQerrorCodeDetail::unableToGetDataBlock       , "unableToGetDataBlock" },
    {CrvDAQerrorCodeDetail::invalidPacket              , "invalidPacket" },
    {CrvDAQerrorCodeDetail::wrongSubsystemID           , "wrongSubsystemID" },
    {CrvDAQerrorCodeDetail::errorUnpackingStatusPacket , "errorUnpackingStatusPacket" },
    {CrvDAQerrorCodeDetail::errorUnpackingCrvHits      , "errorUnpackingCrvHits" },
    {CrvDAQerrorCodeDetail::byteCountMismatch          , "byteCountMismatch" }
  };

  std::map<CrvDAQerrorCodeDetail::enum_type,std::string> const& CrvDAQerrorCodeDetail::names()
  {
    return nam;
  }

}
