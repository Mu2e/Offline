#ifndef CaloConditions_CalSourceCalibConfigConfig_hh
#define CaloConditions_CalSourceCalibConfig_hh
//
// Initialize Calibration of calorimeter from fcl
//
#include <string>
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Sequence.h"

namespace mu2e {

  struct CalSourceCalibConfig {
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    fhicl::Atom<int> verbose{Name("verbose"), Comment("verbosity: 0 or 1")}; 
    fhicl::Atom<bool> useDb{Name("useDb"), Comment("use database or fcl")}; 
    fhicl::Atom<double> roid {Name("roid"), Comment("readout ID")};
    fhicl::Atom<double> EPeak {Name("EPeak"), Comment("peak energy from fit")};
    fhicl::Atom<double> Esource {Name("Esource"), Comment("source energy")};
    fhicl::Atom<double> chisq {Name("chisq"), Comment("chisq of fit")};
  };

}

#endif
