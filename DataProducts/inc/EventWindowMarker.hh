#ifndef DataProducts_EventWindowMarker_hh
#define DataProducts_EventWindowMarker_hh

namespace mu2e {
  struct EventWindowMarker {
    enum SpillType { offspill=0,onspill=1 }; 

    SpillType spillType() const { return _spillType;}
    double eventLength() const { return _eventLength;}
    //
    
    SpillType _spillType; 
    double _eventLength; // in ns, for on spill should be 1675 or 1700
  };
}

#endif
