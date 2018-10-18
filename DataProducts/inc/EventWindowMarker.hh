#ifndef DataProducts_EventWindowMarker_hh
#define DataProducts_EventWindowMarker_hh

namespace mu2e {
  struct EventWindowMarker {
    float timeOffset() const { return _timeOffset;}
    //
    float _timeOffset;
  };
}

#endif
