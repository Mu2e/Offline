#ifndef DataProducts_ExtMonFNALModuleDenseId_hh
#define DataProducts_ExtMonFNALModuleDenseId_hh

// Sequential number of a silicon module in Mu2e ExtMonFNAL detector.
// Zero based.
//
//
// Original author Andrei Gaponenko

#include <ostream>

namespace mu2e {

  class ExtMonFNALModuleDenseId {
  public:

    static const unsigned int NOMODULE = -1u;

    explicit ExtMonFNALModuleDenseId(unsigned int did = NOMODULE) : did_(did) {}

    bool operator==( ExtMonFNALModuleDenseId const& rhs) const{
      return did_ == rhs.did_;
    }

    bool operator!=( ExtMonFNALModuleDenseId const& rhs) const{
      return !(*this == rhs);
    }

    bool operator<( ExtMonFNALModuleDenseId const& rhs) const{
      return did_ < rhs.did_;
    }

    unsigned int number() const { return did_; }

  private:
    unsigned int did_;
  };

  inline std::ostream& operator<<( std::ostream& os, const ExtMonFNALModuleDenseId& id) {
    return os<<id.number();
  }

}
#endif /* DataProducts_ExtMonFNALModuleDenseId_hh */
