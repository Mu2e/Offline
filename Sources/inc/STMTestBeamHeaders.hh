#ifndef RecoDataProducts_STMTestBeamHeaders_hh
#define RecoDataProducts_STMTestBeamHeaders_hh
//
// Data product that represents the digitized signal coming from the STM
//

// C++ includes
#include <iostream>
#include <vector>
#include <array>
#include <Rtypes.h>
#include <ctime>
#include <iomanip>

namespace mu2e {

  namespace STMTestBeam {
    enum TriggerType { kInternal=0, kExternal=1};

    struct TriggerHeader {
      uint16_t deadbeef[2]; // fixed header word (= 0xDEADBEEF)
      uint16_t datasize[2]; // data size
      uint16_t n_slices[2]; // number of slices
      uint16_t trigger_number[2]; // trigger number
      uint16_t trigger_mode[1]; // trigger mode
      uint16_t trigger_time[4]; // trigger time
      uint16_t trigger_offset[2]; // time of first ADC value relative to trigger time
      uint16_t n_dropped_packets[1]; // number of dropped packets
      uint16_t unixtime[4];

      uint32_t getFixedHeader() const {
        return ((uint32_t) deadbeef[1] << 16) | ((uint32_t) deadbeef[0]);
      }

      uint32_t getDataSize() const {
        return ((uint32_t) datasize[1] << 16) | ((uint32_t) datasize[0]);
      }

      uint32_t getNSlices() const {
        return ((uint32_t) n_slices[1] << 16) | ((uint32_t) n_slices[0]);
      }

      uint32_t getTriggerNumber() const {
        return ((uint32_t) trigger_number[1] << 16) | ((uint32_t) trigger_number[0]);
      }

      uint16_t getTriggerMode() const {
        return ((uint16_t) trigger_mode[0]);
      }

      uint64_t getTriggerTime() const {
        return ((uint64_t) trigger_time[3] << 48) | ((uint64_t) trigger_time[2] << 32) | ((uint64_t) trigger_time[1] << 16) | ((uint64_t) trigger_time[0]);
      }

      uint32_t getTriggerOffset() const {
        return ((uint32_t) trigger_offset[1] << 16) | ((uint32_t) trigger_offset[0]);
      }

      uint16_t getNDroppedPackets() const {
        return ((uint16_t) n_dropped_packets[0]);
      }

      uint64_t getUnixTime() const {
        return ((uint64_t) unixtime[3] << 48) | ((uint64_t) unixtime[2] << 32) | ((uint64_t) unixtime[1] << 16) | ((uint64_t) unixtime[0]);
      }

      bool checkFixedHeader() const {
        return (this->getFixedHeader() == (uint32_t) 0xDEADBEEF);
      }


      friend std::ostream& operator<<(std::ostream& os, const TriggerHeader& header) {
        os << "Correct fixed header? " << std::boolalpha << header.checkFixedHeader() << std::endl;
        os << "Data Size: " << header.getDataSize() << " bytes" << std::endl;
        os << "No. of Slices: " << header.getNSlices() << std::endl;
        os << "Trigger #" << header.getTriggerNumber() << std::endl;
        os << "\tMode: " << header.getTriggerMode() << std::endl;
        os << "\tTime: " << header.getTriggerTime() << " ns" << std::endl;
        os << "\tTime Offset: " << header.getTriggerOffset() << " ns" << std::endl;
        os << "\tNo. of Dropped Packets: " << header.getNDroppedPackets() << std::endl;
        using time_point = std::chrono::system_clock::time_point;
        time_point header_timepoint(std::chrono::duration_cast<time_point::duration>(std::chrono::milliseconds(header.getUnixTime())));
        std::time_t header_t = std::chrono::system_clock::to_time_t(header_timepoint);
        os << "Unix Time: " << std::put_time(std::gmtime(&header_t), "%c %Z") << std::endl;
        return os;
      }
    };

    struct SliceHeader {
      uint16_t slice_number[2]; // slice number
      uint16_t slice_size[2]; // slice size
      uint16_t slice_time[4]; // slice time

      uint32_t getSliceNumber() const {
        return ((uint32_t) slice_number[1] << 16) | ((uint32_t) slice_number[0]);
      }

      uint32_t getSliceSize() const {
        return ((uint32_t) slice_size[1] << 16) | ((uint32_t) slice_size[0]);
      }

      uint64_t getSliceTime() const {
        return ((uint64_t) slice_time[3] << 48) | ((uint64_t) slice_time[2] << 32) | ((uint64_t) slice_time[1] << 16) | ((uint64_t) slice_time[0]);
      }

      unsigned long int getNADC() const {
        return (this->getSliceSize())/ sizeof(int16_t);
      }

      friend std::ostream& operator<<(std::ostream& os, const SliceHeader& header) {
        os << "Slice #" << header.getSliceNumber() << std::endl;
        os << "\tSize: " << header.getSliceSize() << " bytes = " << header.getNADC() << " samples" << std::endl;
        os << "\tTime: " << header.getSliceTime() << " ns" << std::endl;

        return os;
      }
    };
  };
}
#endif
