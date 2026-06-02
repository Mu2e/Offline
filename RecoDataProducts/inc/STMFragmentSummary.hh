#ifndef RecoDataProducts_STMFragmentSummary_hh
#define RecoDataProducts_STMFragmentSummary_hh

//
// Data prodcuts that represent the container and inner frag counts
// This is used for the unpacking module
//

// C++ includes
#include <iostream>
#include <vector>
#include <array>
#include <Rtypes.h>

#include "Offline/DataProducts/inc/STMChannel.hh"

namespace mu2e {
  class STMFragmentSummary {

  public:
    STMFragmentSummary(size_t nContainerFrags, size_t nInnerFrags) : _nContainerFrags(nContainerFrags), _nInnerFrags(nInnerFrags) {};

    size_t nContainerFrags() const { return _nContainerFrags; }
    size_t nInnerFrags() const { return _nInnerFrags; }

  private:
    size_t _nContainerFrags{};
    size_t _nInnerFrags{};
  };
  typedef std::vector<STMFragmentSummary> STMFragmentSummaryCollection;
}
#endif
