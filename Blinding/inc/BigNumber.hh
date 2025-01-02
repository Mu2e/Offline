// Ed Callaghan
// Interface for a large number; stored as a string, to be manipulated with gmp
// August 2024

#ifndef BigNumber_hh
#define BigNumber_hh

#include <string>

namespace mu2e{
  class BigNumber{
    public:
      BigNumber() = default;
      BigNumber(std::string);
      const char* Buffer() const;
      const std::string& String() const;
      bool IsZero() const;
      bool operator== (const BigNumber&) const;
    protected:
      std::string _rep;
    private:
      /**/
  };
} // namespace mu2e

#endif
