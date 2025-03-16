#include <type_traits>
#include <utility>

#include "Offline/RecoDataProducts/inc/DAQerror.hh"

namespace mu2e
{
  std::string const& DAQerrorCodeDetail::typeName()
  {
    static const std::string type("DAQerrorCode");
    return type;
  }

  static const std::map<DAQerrorCodeDetail::enum_type,std::string> nam
  {
    {DAQerrorCodeDetail::unknown                    , "unknown" },
    {DAQerrorCodeDetail::byteCountMismatch          , "byteCountMismatch" }
  };

  std::map<DAQerrorCodeDetail::enum_type,std::string> const& DAQerrorCodeDetail::names()
  {
    return nam;
  }

}
