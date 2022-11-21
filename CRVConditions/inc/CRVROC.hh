#ifndef CRVConditions_CRVROC_hh
#define CRVConditions_CRVROC_hh

//
// a little holder for the triplet of ROC numbering
//

#include <cstddef>

namespace mu2e {

class CRVROC {
 public:
  CRVROC(std::size_t ROC, std::size_t FEB, std::size_t FEBchannel) :
      _ROC(ROC), _FEB(FEB), _FEBchannel(FEBchannel) {}
  std::size_t ROC() const { return _ROC; }
  std::size_t FEB() const { return _FEB; }
  std::size_t FEBchannel() const { return _FEBchannel; }
  std::size_t _ROC;
  std::size_t _FEB;
  std::size_t _FEBchannel;
};

}  // namespace mu2e

#endif
