//
// Trivial typedefs to define generic 'RecoQual' quality value from reco algorithms
//  All blame goes to: Dave Brown (LBNL) 3 June 2019
//
#ifndef RecoDataProducts_RecoQual_hh
#define RecoDataProducts_RecoQual_hh
#include <Rtypes.h>
#include <vector>
#include "GeneralUtilities/inc/MVAStatus.hh"
namespace mu2e {
  struct RecoQual {
    RecoQual() : _status(MVAStatus::unset), _value(-1.0) {}
    RecoQual(MVAStatus status, Float_t value ) : _status(status), _value(value) {}
    MVAStatus _status;
    Float_t _value;
  };
  typedef std::vector<RecoQual> RecoQualCollection;
}
#endif
