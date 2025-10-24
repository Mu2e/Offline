#ifndef DataProducts_CRVId_hh
#define DataProducts_CRVId_hh
//
// Currenlty, CRV bar count is determined in the process
// of building the geometry.  This constant
// which should equal the generated count.
// It is used in the conditions where the count can't be dynamic
// TODO  When the CRV geometry is final, this number will have to change

#include <cstddef>

namespace mu2e {

  class CRVId{
  public:

    // taken from dynamic geometry
    constexpr static std::size_t nBars = 5504;
    constexpr static std::size_t nChanPerBar = 4;
    constexpr static std::size_t nSidesPerBar = 2;
    // some bars have only 2 SiPMs so this count is sparse
    constexpr static std::size_t nChannels = nBars*nChanPerBar;
    constexpr static std::size_t nLayers = 4;
    // not all the implied possible channels below are active
    // these are the dimensions of sparse containers
    constexpr static std::size_t nChanPerFEB = 64;
    constexpr static std::size_t nFEBPerROC = 25;
    constexpr static std::size_t nROC = 18;
    constexpr static std::size_t nROCPerDTC = 6;

    constexpr static std::size_t nChanPerFPGA = 16;
    constexpr static std::size_t nFPGAPerFEB = 4;
  };
}
#endif /* DataProducts_CRVId_hh */
