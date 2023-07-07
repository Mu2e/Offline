#ifndef DataProducts_STMTestBeamHeaders_hh
#define DataProducts_STMTestBeamHeaders_hh
//
// Header definitions for STM test beam data
//

// C++ includes
#include <iostream>
#include <vector>
#include <array>
#include <ctime>
#include <iomanip>

namespace mu2e {

  namespace STMTestBeam {
    struct TriggerHeader {
      uint16_t deadbeef[2]; // fixed header word (= 0xDEADBEEF)
      uint16_t datasize[2]; // data size
      uint16_t n_slices[2]; // number of slices
      uint16_t trigger_number[2]; // trigger number
      uint16_t trigger_mode[1]; // trigger mode
      uint16_t trigger_time[4]; // trigger time
      uint16_t trigger_offset[2]; // time of first ADC value relative to trigger time
      uint16_t n_dropped_packets[1]; // number of dropped packets

      uint32_t getFixedHeader() const {
        return ( uint32_t(deadbeef[1]) << 16 | uint32_t(deadbeef[0]) );
      }

      uint32_t getDataSize() const {
        return ( uint32_t(datasize[1]) << 16 | uint32_t(datasize[0]) );
      }

      uint32_t getNSlices() const {
        return ( uint32_t(n_slices[1]) << 16 | uint32_t(n_slices[0]) );
      }

      uint32_t getTriggerNumber() const {
        return ( uint32_t(trigger_number[1]) << 16 | uint32_t(trigger_number[0]) );
      }

      uint16_t getTriggerMode() const {
        return uint16_t(trigger_mode[0]);
      }

      uint64_t getTriggerTime() const {
        return ( uint64_t(trigger_time[3]) << 48 | uint64_t(trigger_time[2]) << 32 | uint64_t(trigger_time[1]) << 16 | uint64_t(trigger_time[0]) );
      }

      uint32_t getTriggerOffset() const {
        return ( uint32_t(trigger_offset[1]) << 16 | uint32_t(trigger_offset[0]) );
      }

      uint16_t getNDroppedPackets() const {
        return uint16_t(n_dropped_packets[0]);
      }

      bool checkFixedHeader() const {
        return (this->getFixedHeader() == uint32_t(0xDEADBEEF));
      }
    };

    struct SliceHeader {
      uint16_t slice_number[2]; // slice number
      uint16_t slice_size[2]; // slice size
      uint16_t slice_time[4]; // slice time

      uint32_t getSliceNumber() const {
        return ( uint32_t(slice_number[1]) << 16 | uint32_t(slice_number[0]) );
      }

      uint32_t getSliceSize() const {
        return ( uint32_t(slice_size[1]) << 16 | uint32_t(slice_size[0]) );
      }

      uint64_t getSliceTime() const {
        return ( uint64_t(slice_time[3]) << 48 | uint64_t(slice_time[2]) << 32 | uint64_t(slice_time[1]) << 16 | uint64_t(slice_time[0]) );
      }

      unsigned long int getNADC() const {
        return (this->getSliceSize())/ sizeof(int16_t);
      }
    };
  }

  std::ostream& operator<<(std::ostream& os, const STMTestBeam::TriggerHeader& header);
  std::ostream& operator<<(std::ostream& os, const STMTestBeam::SliceHeader& header);
}
#endif
