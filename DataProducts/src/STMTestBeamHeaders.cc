//
// Stream operators for STMTestBeamHeaders
//

// C++ includes
#include "Offline/DataProducts/inc/STMTestBeamHeaders.hh"

namespace mu2e {

  std::ostream& operator<<(std::ostream& os, const STMTestBeam::TriggerHeader& header) {
    os << "Correct fixed header? " << std::boolalpha << header.checkFixedHeader() << std::endl;
    os << "Data Size: " << header.getDataSize() << " bytes" << std::endl;
    os << "No. of Slices: " << header.getNSlices() << std::endl;
    os << "Trigger #" << header.getTriggerNumber() << std::endl;
    os << "\tMode: " << header.getTriggerMode() << std::endl;
    os << "\tTime: " << header.getTriggerTime() << " ct" << std::endl;
    os << "\tTime Offset: " << header.getTriggerOffset() << " ct" << std::endl;
    os << "\tNo. of Dropped Packets: " << header.getNDroppedPackets() << std::endl;
    return os;
  }

  std::ostream& operator<<(std::ostream& os, const STMTestBeam::SliceHeader& header) {
    os << "Slice #" << header.getSliceNumber() << std::endl;
    os << "\tSize: " << header.getSliceSize() << " bytes = " << header.getNADC() << " samples" << std::endl;
    os << "\tTime: " << header.getSliceTime() << " ct" << std::endl;

    return os;
  }
}
