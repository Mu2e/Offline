#ifndef DataProducts_CaloRawId_hh
#define DataProducts_CaloRawId_hh
//
//
// Online or raw data identifier of one calorimeter SiPM channel
// Offline numbering is in CaloId
//
#include <iosfwd>
#include <string>
#include <math.h>

namespace mu2e {

  class CaloRawId{

  public:
    constexpr static uint16_t _nDIRAC       = 136;
    constexpr static uint16_t _nChPerDIRAC  = 20;
    constexpr static uint16_t _nChannel     = _nChPerDIRAC*_nDIRAC;

    explicit CaloRawId(uint16_t id) { _id = id; }

    uint16_t channel() { return _id; }
    uint16_t dirac() { return _id/20; }
    uint16_t ROCchannel() { return _id%20; }

    bool isValid() { return _id < _nChannel; }

    uint16_t _id;
  };

}
#endif /* DataProducts_CaloRawId_hh */
