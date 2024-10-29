#ifndef CaloConditions_CalCalibConfig_hh
#define CaloConditions_CalCalibConfig_hh
//
// Initialize calibration info from fcl
// author: S. Middleton 2022
#include <string>
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Sequence.h"

namespace mu2e {

  struct CalCalibConfig {
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    fhicl::Atom<int> verbose{Name("verbose"), Comment("verbosity: 0 or 1")};
    fhicl::Atom<bool> useDb{Name("useDb"), Comment("use database or fcl")};
    //fhicl::Atom<uint16_t> roid {Name("roid"), Comment("unique offline readout ID")};
    fhicl::Atom<float> ADC2MeV {Name("ADC2MeV"), Comment("constant per SiPM")};
    fhicl::Atom<float> timeoffset {Name("timeoffset"), Comment("constant per SiPM")};
  };

}

#endif
