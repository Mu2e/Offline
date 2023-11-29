#ifndef DAQConditions_EventTimingConfig_hh
#define DAQConditions_EventTimingConfig_hh
//
// Initialize EventTiming from fcl
//
#include <string>
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

namespace mu2e {

  struct EventTimingConfig {
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;

    fhicl::Atom<int> verbose{
      Name("verbose"), Comment("verbosity: 0 or 1")};
    fhicl::Atom<bool> useDb{
      Name("useDb"), Comment("use database or fcl")};

    fhicl::Atom<double> timeFromProtonsToDRMarker{
      Name("TimeFromProtonsToDRMarker"), Comment("Time shift in ns of DR marker wrt to proton peak. Positive means marker arrives after protons")};
    fhicl::Atom<unsigned> offSpillLength{
      Name("OffSpillEventLength"), Comment("Length of off spill events in clock counts" )};

  };

}

#endif
