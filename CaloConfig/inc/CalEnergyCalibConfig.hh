#ifndef CaloConditions_CalEnergyCalibConfig_hh
#define CaloConditions_CalEnergyCalibConfig_hh
//
// Initialize Calibration of calorimeter from fcl
//
#include <string>
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Sequence.h"

namespace mu2e {

  struct CalEnergyCalibConfig {
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    fhicl::Atom<int> verbose{Name("verbose"), Comment("verbosity: 0 or 1")}; 
    fhicl::Atom<bool> useDb{Name("useDb"), Comment("use database or fcl")}; 
    fhicl::Atom<double> roid {Name("roid"), Comment("readout ID")};
};

}

#endif
