#ifndef DataProducts_CaloId_hh
#define DataProducts_CaloId_hh
//
// Identifier of one straw in a tracker.
// Original author Rob Kutschke
// Re-implemented as integer bitfields by Dave Brown (LBNL)
//
#include <iosfwd>
#include <string>
#include <math.h>

namespace mu2e {

  class CaloId{

    // define the bit field shifts and masks
  public:
    constexpr static uint16_t _nDIRAC          = 136;
    constexpr static uint16_t _nChPerDIRAC     = 20;
    constexpr static uint16_t _nCrystalPerDisk = 674;
    constexpr static uint16_t _nSiPMPerCrystal = 2;
    constexpr static uint16_t _nDisk           = 2;
    constexpr static uint16_t _nCrystalCAPHRI  = 4;
    constexpr static uint16_t _nPINDiodPerDisk = 8;
    constexpr static uint16_t _nTotChannel     = _nChPerDIRAC*_nDIRAC;
    constexpr static uint16_t _nCrystalChannel = _nCrystalPerDisk*_nSiPMPerCrystal; //total number of readout channels associated to the CsI and CAPHRI crystals

    CaloId(){}

  };

}
#endif /* DataProducts_CaloId_hh */
