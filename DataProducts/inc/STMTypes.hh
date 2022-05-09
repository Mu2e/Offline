//
// Types used by STM data products
//
#ifndef DataProducts_STMTypes_hh_
#define DataProducts_STMTypes_hh_

namespace mu2e {
  enum STMTriggerMode { kInternal=0, kExternal=1 };
  enum STMChannel { kHPGe=0, kLaBr=1 };
  enum STMDataType { kUnsuppressed=0, kZeroSuppressed=1, kMWD=2};

  struct STMTrigType {
  public:
    STMTrigType(uint16_t mode, uint16_t channel, uint16_t type) {
      _data = (mode << 12) | (channel << 10) | type;
    }

    STMTriggerMode mode() const { return STMTriggerMode((_data & 0xF000) >> 12); }
    STMChannel channel() const { return STMChannel((_data & 0xC00) >> 10); }
    STMDataType type() const { return STMDataType((_data & 0x3FF)); }
  private:
    uint16_t _data;
  };
}

#endif
