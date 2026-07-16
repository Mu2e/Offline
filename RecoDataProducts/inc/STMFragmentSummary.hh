#ifndef RecoDataProducts_STMFragmentSummary_hh
#define RecoDataProducts_STMFragmentSummary_hh
//
// Data products containing STM  container and inner fragment counts
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
    STMFragmentSummary() : _nContainerFrags(0), _nInnerFrags(0),
                           _nGoodRawFrags(0),_nGoodZSFrags(0),_nGoodPHFrags(0),
                           _nZeroRawFrags(0),_nZeroZSFrags(0),_nZeroPHFrags(0),
                           _nEmptyRawFrags(0),_nEmptyZSFrags(0),_nEmptyPHFrags(0){};

    STMFragmentSummary(size_t nContainerFrags, size_t nInnerFrags,
                       size_t nGoodRawFrags,size_t nGoodZSFrags, size_t nGoodPHFrags,
                       size_t nZeroRawFrags, size_t nZeroZSFrags, size_t nZeroPHFrags,
                       size_t nEmptyRawFrags, size_t nEmptyZSFrags, size_t nEmptyPHFrags) :
      _nContainerFrags(nContainerFrags), _nInnerFrags(nInnerFrags),
      _nGoodRawFrags(nGoodRawFrags),_nGoodZSFrags(nGoodZSFrags),_nGoodPHFrags(nGoodPHFrags),
      _nZeroRawFrags(nZeroRawFrags),_nZeroZSFrags(nZeroZSFrags),_nZeroPHFrags(nZeroPHFrags),
      _nEmptyRawFrags(nEmptyRawFrags),_nEmptyZSFrags(nEmptyZSFrags),_nEmptyPHFrags(nEmptyPHFrags) {};


    size_t nContainerFrags() const { return _nContainerFrags; }
    size_t nInnerFrags() const { return _nInnerFrags; }

    size_t nGoodRawFrags() const { return _nGoodRawFrags; }
    size_t nGoodZSFrags() const { return _nGoodZSFrags; }
    size_t nGoodPHFrags() const { return _nGoodPHFrags; }

    size_t nZeroRawFrags() const { return _nZeroRawFrags; }
    size_t nZeroZSFrags() const { return _nZeroZSFrags; }
    size_t nZeroPHFrags() const { return _nZeroPHFrags; }

    size_t nEmptyRawFrags() const { return _nEmptyRawFrags; }
    size_t nEmptyZSFrags() const { return _nEmptyZSFrags; }
    size_t nEmptyPHFrags() const { return _nEmptyPHFrags; }


  private:
    size_t _nContainerFrags{};
    size_t _nInnerFrags{};

    size_t _nGoodRawFrags{};
    size_t _nGoodZSFrags{};
    size_t _nGoodPHFrags{};

    size_t _nZeroRawFrags{};
    size_t _nZeroZSFrags{};
    size_t _nZeroPHFrags{};

    size_t _nEmptyRawFrags{};
    size_t _nEmptyZSFrags{};
    size_t _nEmptyPHFrags{};


  };
  typedef std::vector<STMFragmentSummary> STMFragmentSummaryCollection;
}
#endif
