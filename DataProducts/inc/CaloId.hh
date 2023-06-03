#ifndef DataProducts_CaloId_hh
#define DataProducts_CaloId_hh
//
// Offline identifier of one calorimeter SiPM channel
// Online numbering is in CaloRawId
//
#include <array>
#include <algorithm>

namespace mu2e {

  class CaloId{

    // define the bit field shifts and masks
  public:
    constexpr static uint16_t _nCrystalPerDisk = 674;
    constexpr static uint16_t _nSiPMPerCrystal = 2;
    constexpr static uint16_t _nDisk           = 2;
    // PIN diodes are each one SiPM channel
    constexpr static uint16_t _nPINDiodPerDisk = 8;
    constexpr static uint16_t _nCrystal        = _nCrystalPerDisk*_nDisk;
    constexpr static uint16_t _nCrystalChannel = _nCrystalPerDisk*_nSiPMPerCrystal;
    constexpr static uint16_t _nCaphriCrystal  = 4;
    // crystal numbers, not SiPM channels. Only in disk 0.
    constexpr static std::array<uint16_t,_nCaphriCrystal> _caphriId = {582,609,610,637};
    constexpr static uint16_t _nChannel        = _nCrystalChannel + _nPINDiodPerDisk*_nDisk;

    explicit CaloId(uint16_t id) { _id = id; }

    uint16_t channel() { return _id; }
    uint16_t SiPM01() { return _id%_nSiPMPerCrystal; }
    uint16_t crystal() { return _id/_nSiPMPerCrystal; }
    uint16_t disk();
    bool isValid() { return _id < _nChannel; }
    bool isCrystal() { return _id < _nCrystal; }
    bool isCaphri();
    bool isPINDiode() { return _id >= _nCrystal; }

    uint16_t _id;
  };

};
#endif /* DataProducts_CaloId_hh */
