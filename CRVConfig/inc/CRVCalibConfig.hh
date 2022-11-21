#ifndef STMConfig_CRVCalibConfig_hh
#define STMConfig_CRVCalibConfig_hh
//
// Fcl stanza for CRV pedestal and energy calibration prodition
//
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Tuple.h"
#include <string>

namespace mu2e {

struct CRVCalibConfig {
  using Name = fhicl::Name;
  using Comment = fhicl::Comment;

  fhicl::Atom<int> verbose{Name("verbose"), Comment("verbosity: 0 or 1")};
  fhicl::Atom<bool> useDb{Name("useDb"), Comment("use database or fcl")};

  fhicl::Atom<float> pedestal{Name("pedestal"), Comment("universal pedestal")};
  fhicl::Atom<float> height{Name("height"), Comment("universal height")};
  fhicl::Atom<float> area{Name("area"), Comment("universal area")};
  fhicl::Atom<float> timeOffset{Name("timeOffset"),
                                Comment("universal timeOffset")};
};

}  // namespace mu2e

#endif
