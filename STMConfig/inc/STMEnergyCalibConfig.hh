#ifndef STMConfig_STMEnergyCalibConfig_hh
#define STMConfig_STMEnergyCalibConfig_hh
//
// Fcl stanza for STM energy prodition
//
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Tuple.h"
#include <string>

namespace mu2e {

struct STMEnergyCalibConfig {
  using Name = fhicl::Name;
  using Comment = fhicl::Comment;

  fhicl::Atom<int> verbose{Name("verbose"), Comment("verbosity: 0 or 1")};
  fhicl::Atom<bool> useDb{Name("useDb"), Comment("use database or fcl")};
  fhicl::Sequence<fhicl::Tuple<std::string, float>, 2> pedestals{
      Name("pedestals"),
      Comment("sequence of tuples of channel_name and value")};
};

}  // namespace mu2e

#endif
