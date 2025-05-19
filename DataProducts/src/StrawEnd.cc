#include "Offline/DataProducts/inc/StrawEnd.hh"
namespace mu2e {
  std::ostream& operator << (std::ostream& os, mu2e::StrawEnd const& end){
    switch ( end._end ) {
      case mu2e::StrawEnd::cal:
        os << "Cal";
        break;
      case mu2e::StrawEnd::hv:
        os << "HV";
        break;
      default:
        os << "Unknown";
        break;
    }
    return os;
  }
}
