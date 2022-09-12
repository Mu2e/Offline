//
// Types used by STM data products
//
#ifndef DataProducts_STMTypes_hh_
#define DataProducts_STMTypes_hh_

#include "Offline/DataProducts/inc/STMChannel.hh"

namespace mu2e {
  enum STMTriggerMode { kInternal=0, kExternal=1 };
  enum STMDataType { kUnsuppressed=0, kZeroSuppressed=1, kMWD=2, kPQ=3};

  struct STMTrigType {
  public:
    STMTrigType(uint16_t mode, uint16_t channel, uint16_t type) {
      _data = (mode) | (channel << 10) | type;
    }
    uint16_t data() const { return _data; }
    STMTriggerMode mode() const { return STMTriggerMode((_data & 0xF000) >> 12); }
    STMChannel channel() const {
      STMChannel::enum_type ch = static_cast<STMChannel::enum_type>((_data & 0xC00) >> 10);
      return STMChannel(ch);
    }
    STMDataType type() const { return STMDataType((_data & 0x3FF)); }
  private:
    uint16_t _data;
  };

  enum STMDigiFlag : uint16_t { kUnknown = 0, kOK = 1 };
}

#endif
