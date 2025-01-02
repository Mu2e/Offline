// Ed Callaghan
// Interface for a large number; stored as a string, to be manipulated with gmp
// August 2024

#include "Offline/Blinding/inc/BigNumber.hh"

namespace mu2e{
  BigNumber::BigNumber(std::string rep): _rep(rep){
    /**/
  }

  const char* BigNumber::Buffer() const{
    auto rv = _rep.c_str();
    return rv;
  }

  const std::string& BigNumber::String() const{
    auto& rv = _rep;
    return rv;
  }

  bool BigNumber::IsZero() const{
    auto rv = (_rep == "0");
    return rv;
  }

  bool BigNumber::operator== (const BigNumber& rhs) const{
    auto rv = (_rep == rhs.String());
    return rv;
  }
} // namespace mu2e
