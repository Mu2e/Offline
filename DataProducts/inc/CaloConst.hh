#ifndef DataProducts_CaloConst_hh
#define DataProducts_CaloConst_hh
//
// Constants describing the calorimeter
// used in geometry, channel number ID classes, and database access
//
#include <array>
#include <cstdint>

namespace mu2e {

  class CaloConst{

    // define the bit field shifts and masks
  public:
    constexpr static uint16_t _nCrystalPerDisk = 674;
    constexpr static uint16_t _nSiPMPerCrystal = 2;
    constexpr static uint16_t _nDisk           = 2;
    // PIN diodes are each one SiPM channel
    constexpr static uint16_t _nPINDiodPerDisk = 8;
    constexpr static uint16_t _nPINDiodeLaserBox = 4;
    constexpr static uint16_t _nChannelSpares = 4;
    constexpr static uint16_t _nCrystal        = _nCrystalPerDisk*_nDisk;
    constexpr static uint16_t _nCrystalChannel = _nCrystal*_nSiPMPerCrystal;
    constexpr static uint16_t _nCaphriCrystal  = 4;
    // crystal numbers, not SiPM channels. Only in disk 0.
    constexpr static std::array<uint16_t,_nCaphriCrystal> _caphriId = {582,609,610,637};

    constexpr static uint16_t _nChannel        = _nCrystalChannel + _nPINDiodPerDisk*_nDisk + _nPINDiodeLaserBox;
    constexpr static uint16_t _nChannelDB      = _nChannel + _nChannelSpares;

    constexpr static uint16_t _nDIRAC       = 161;
    constexpr static uint16_t _nChPerDIRAC  = 20;
    constexpr static uint16_t _nRawChannel  = _nChPerDIRAC*_nDIRAC;

    constexpr static uint16_t _invalid      = 9999;

    using CaloSiPMId_type = std::uint16_t;

    enum SiPMn {SiPM0=0,SiPM1=1};
    enum detType {CsI=0,CAPHRI=1,PINDiode=2,Invalid=3};

  };

}
#endif /* DataProducts_CaloConst_hh */
