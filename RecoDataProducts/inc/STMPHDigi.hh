#ifndef RecoDataProducts_STMPHDigi_hh
#define RecoDataProducts_STMPHDigi_hh
//
// Data product that represents an uncalibrated hit in the STM
//

// C++ includes
#include <iostream>
#include <vector>
#include <array>
#include <Rtypes.h>

namespace mu2e {

  class STMPHDigi {
  public:
    STMPHDigi() : _time(0), _energy(0){};

    STMPHDigi(uint32_t time, int16_t energy) : _time(time), _energy(energy) {};


    uint32_t time() const { return _time; }
    int16_t energy() const { return _energy; }

  private:
    uint32_t _time; // uncalibrated time of hit [ct]
    int16_t _energy; // uncalibrated energy [ADC units]
  };

  typedef std::vector<STMPHDigi> STMPHDigiCollection;

  bool lessByTime(const STMPHDigi& a, const STMPHDigi& b) {
    if (a.time() < b.time()) { return true; }
    else { return false; }
  }

  bool lessByEnergy(const STMPHDigi& a, const STMPHDigi& b) {
    if (a.energy() < b.energy()) { return true; }
    else { return false; }
  }
}
#endif
