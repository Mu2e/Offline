#ifndef RecoDataProducts_STMTestBeamPacketDefns_hh
#define RecoDataProducts_STMTestBeamPacketDefns_hh
//
// Data product that represents the digitized signal coming from the STM
//

// C++ includes
#include <iostream>
#include <vector>
#include <array>
#include <Rtypes.h>

namespace mu2e {

  enum STMTestBeamTriggerMode { kExternal=0, kInternal=1};

  struct STMTestBeamTriggerHeader {
    uint16_t macropulseNumber[2]; // macropulse number
    uint16_t macropulseTime[4]; // time of the macropulse
    uint16_t triggerMode[1]; // trigger mode
    uint16_t microCount[2]; // ???
    uint16_t channel[1]; // channel number
    uint16_t offset[1]; // ???
    uint16_t nSlices[1]; // number of slices per packet
    uint16_t unused[4]; // unused

    uint32_t getMacropulseNumber() const {
      return (uint32_t)macropulseNumber[0] << 16 | (uint32_t)macropulseNumber[1];
    }

    uint64_t getMacropulseTime() const {
      return (uint64_t)macropulseTime[0] << 48 | (uint64_t)macropulseTime[1] << 32 | (uint64_t)macropulseTime[2] << 16 | (uint64_t)macropulseTime[3];
    }

    uint16_t getTriggerMode() const {
      return triggerMode[0];
    }

    uint32_t getMicroCount() const {
      return (uint32_t)microCount[0] << 16 | (uint32_t)microCount[1];
    }

    uint16_t getChannel() const {
      return channel[0];
    }

    uint16_t getOffset() const {
      return offset[0];
    }

    uint16_t getNSlices() const {
      return nSlices[0];
    }


    friend std::ostream& operator<<(std::ostream& os, const STMTestBeamTriggerHeader& triggerHeader) {
      os << "Macropulse #" << triggerHeader.getMacropulseNumber() << std::endl;
      os << "\tTime: " << triggerHeader.getMacropulseTime() << " ns" << std::endl;
      os << "\tTrigger Mode: " << triggerHeader.getTriggerMode() << std::endl;
      os << "\tMicro Count: " << triggerHeader.getMicroCount() << std::endl;
      os << "\tChannel: " << triggerHeader.getChannel() << std::endl;
      os << "\tOffset: " << triggerHeader.getOffset() << std::endl;
      os << "\tSlices per Packet: " << triggerHeader.getNSlices() << std::endl;
      return os;
    }
  };

  struct STMTestBeamSliceHeader {
    uint16_t sliceNumber[2]; // slice number
    uint16_t sliceTime[4]; // time of first sample in slice
    uint16_t nSamples[2]; // number of samples in slice

    uint32_t getSliceNumber() const {
      return (uint32_t)sliceNumber[0] << 16 | (uint32_t)sliceNumber[1];
    }

    uint64_t getSliceTime() const {
      return (uint64_t)sliceTime[0] << 48 | (uint64_t)sliceTime[1] << 32 | (uint64_t)sliceTime[2] << 16 | (uint64_t)sliceTime[3];
    }

    uint32_t getNSamples() const {
      return (uint32_t)nSamples[0] << 16 | (uint32_t)nSamples[1];
    }


    friend std::ostream& operator<<(std::ostream& os, const STMTestBeamSliceHeader& sliceHeader) {
      os << "Slice #" << sliceHeader.getSliceNumber() << std::endl;
      os << "\tTime: " << sliceHeader.getSliceTime() << " ns" << std::endl;
      os << "\tNumber of Samples: " << sliceHeader.getNSamples() << std::endl;
      return os;
    }
  };


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
