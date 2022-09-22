#ifndef DataProducts_STMTestBeamEventInfo_hh
#define DataProducts_STMTestBeamEventInfo_hh

namespace mu2e {
  class STMTestBeamEventInfo {
  public:
    enum TrigType { unknown=-1, internal=0, external=1 };
    STMTestBeamEventInfo() : _trigType(TrigType::unknown), _triggerTime(0) { }

    void trigType(uint16_t type) { _trigType = TrigType(type);}
    void triggerTime(uint64_t time) { _triggerTime = time;}

    TrigType trigType() const { return _trigType;}
    uint64_t triggerTime() const { return _triggerTime;}

  private:
    TrigType _trigType;
    uint64_t _triggerTime;
  };
}

#endif
