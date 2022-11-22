#ifndef DataProducts_CRVId_hh
#define DataProducts_CRVId_hh
//
// Currenlty, CRV bar count is determined in the process
// of building the geometry.  This is a constant
// which should equal the generated count.
// It is used in the conditions where the count shouldn't be dynamic
// TODO
// When the CRV geometry is final, this number will have to change

#include <cstddef>

namespace mu2e {

  class CRVId{
  public:

    // taken from dynamic geometry
    constexpr static std::size_t nBars = 5504;
    // some bars have only 2 SiPMs so this is a little sparse
    constexpr static std::size_t nChannels = nBars*4;
    constexpr static std::size_t nChanPerBar = 4;
    constexpr static std::size_t nChanPerFEB = 64;
    constexpr static std::size_t nFEBPerROC = 24;
    constexpr static std::size_t nROC = nChannels/nFEBPerROC/nChanPerFEB + 1;

  };
}
#endif /* DataProducts_CRVId_hh */
