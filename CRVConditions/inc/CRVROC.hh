#ifndef CRVConditions_CRVROC_hh
#define CRVConditions_CRVROC_hh

//
// a little holder for the triplet of ROC numbering
//

#include <cstdint>

namespace mu2e {

class CRVROC {
 public:
  CRVROC(std::uint16_t ROC, std::uint16_t FEB, std::uint16_t FEBchannel) :
      _ROC(ROC), _FEB(FEB), _FEBchannel(FEBchannel) {}
  std::uint16_t ROC() const { return _ROC; }
  std::uint16_t FEB() const { return _FEB; }
  std::uint16_t FEBchannel() const { return _FEBchannel; }
  std::uint16_t _ROC;
  std::uint16_t _FEB;
  std::uint16_t _FEBchannel;
};

}  // namespace mu2e

#endif
