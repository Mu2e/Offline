#ifndef DataProducts_CaloId_hh
#define DataProducts_CaloId_hh
//
// Identifier of one crystal in calorimeter
// Original author Rob Kutschke
// Re-implemented as integer bitfields by Dave Brown (LBNL)
// Edits by S. Middleton (2022)
#include <iosfwd>
#include <string>
#include <math.h>

namespace mu2e {

  class CaloId{

    // define the bit field shifts and masks
  public:
    constexpr static uint16_t _channelmsk = 0x7F;
    constexpr static uint16_t _crystalmsk = 0x7F;
    constexpr static uint16_t _diskmsk = 0x7F;
    constexpr static uint16_t _crystalchannelmsk = 0x7F;

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
    static bool validChannel(uint16_t ichannel) { return ichannel < _nChPerDIRAC;}
    static bool validCrystal (uint16_t icrystal) { return icrystal < _nCrystalPerDisk;}
    static bool validDisk(uint16_t idisk) { return idisk < _nDisk;}
    static bool validCrystalChannel(uint16_t icrych) { return icrych < _nCrystalChannel;}
    bool valid() const { return validChannel(getChannel()) && validCrystal(getCrystal()) && validDisk(getDisk()) && validCrystalChannel(getCrystalChannel());}

    CaloId( uint16_t channel,
      uint16_t crystal,
      uint16_t disk,
      uint16_t crychanel);

    explicit CaloId(){}

    uint16_t getChannel() const{
      return (_sid & _channelmsk);
    }

    uint16_t getCrystal() const{
      return (_sid & _crystalmsk);
    }

    uint16_t getDisk() const{
      return (_sid & _diskmsk);
    }
      
    uint16_t getCrystalChannel() const{
      return (_sid & _crystalchannelmsk);
    }

   private:
      //  data member is a short
      uint16_t _sid;
      // fill fields
      void setChannel(uint16_t ich);
      void setCrystal(uint16_t icry);
      void setDisk(uint16_t id);
      void setCrystalChannel(uint16_t icc);


  };

}
#endif /* DataProducts_CaloId_hh */
