#ifndef DataProducts_CaloId_hh
#define DataProducts_CaloId_hh
//
// Identifier of one crystal in calorimeter
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
    // check that values make sense:
    static bool validChannel(uint16_t ichannel) { return ichannel < _nChPerDIRAC; }
    static bool validCrystal (uint16_t icrystal) { return icrystal < _nCrystalPerDisk; }
    static bool validDisk(uint16_t idisk) { return idisk < _nDisk; }
    static bool validCrystalChannel(uint16_t icrych) { return icrych < _nCrystalChannel; }
      
    CaloId(){}
    
  };

}
#endif /* DataProducts_CaloId_hh */
