#ifndef RecoDataProducts_STMTestBeamBinaryPacket_hh
#define RecoDataProducts_STMTestBeamBinaryPacket_hh
//
// Data product that represents the digitized signal coming from the STM
//

// C++ includes
#include <iostream>
#include <vector>
#include <array>
#include <Rtypes.h>

namespace mu2e {

  struct STMTestBeamBinaryPacket {
    uint16_t trigNum;
    uint16_t trigType;
    uint32_t trigTime;
    uint16_t trigTimeOffset;
    uint16_t baselineMean;
    uint16_t baselineRMS;
    uint16_t nDrop;
    uint16_t ADC0;
    uint16_t ADC1;
    uint16_t ADC2;
    uint16_t ADC3;

    friend std::ostream& operator<<(std::ostream& os, const STMTestBeamBinaryPacket& test) {
      os << "Trigger #" << test.trigNum << std::endl;
      os << "\tTrigger Type: " << test.trigType << std::endl;
      os << "\tTrigger Time: " << test.trigTime << std::endl;
      os << "\tTrigger Time Offset: " << test.trigTimeOffset << std::endl;
      os << "\tBaseline Mean: " << test.baselineMean << std::endl;
      os << "\tBaseline RMS: " << test.baselineRMS << std::endl;
      os << "\tNo. of Dropped Packets: " << test.nDrop << std::endl;
      os << "\tADCs: " << test.ADC0 << " " << test.ADC1 << " " << test.ADC2 << " " << test.ADC3 << std::endl;
      return os;
    }
  };
}
#endif
