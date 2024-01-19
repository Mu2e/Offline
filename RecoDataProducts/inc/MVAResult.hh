//
// A data product that stores the result of any MVA/ML/AI algorithm
//
#ifndef RecoDataProducts_MVAResult_hh
#define RecoDataProducts_MVAResult_hh
#include <Rtypes.h>
#include <vector>

namespace mu2e {
  struct MVAResult {
    MVAResult() : _value(-1.0) {}
    MVAResult(Float_t value ) : _value(value) {}
    Float_t _value;
  };
  typedef std::vector<MVAResult> MVAResultCollection;
}
#endif
