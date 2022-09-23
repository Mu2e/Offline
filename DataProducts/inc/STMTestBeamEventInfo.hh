#ifndef DataProducts_STMTestBeamEventInfo_hh
#define DataProducts_STMTestBeamEventInfo_hh

namespace mu2e {
  class STMTestBeamEventInfo {
  public:
    enum TriggerType { unknown=-1, internal=0, external=1 };
    STMTestBeamEventInfo() { }

    void triggerType(uint16_t type) { _triggerType = TriggerType(type);}
    void triggerTime(uint64_t time) { _triggerTime = time;}

    TriggerType triggerType() const { return _triggerType;}
    uint64_t triggerTime() const { return _triggerTime;}

  private:
    TriggerType _triggerType = TriggerType::unknown;
    uint64_t _triggerTime = 0;
  };
}

#endif
