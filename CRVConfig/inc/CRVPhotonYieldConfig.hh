#ifndef STMConfig_CRVPhotonYieldConfig_hh
#define STMConfig_CRVPhotonYieldConfig_hh
//
// Fcl stanza for CRV PhotonYield (scintillation yield spread of CRV channels) prodition
//
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Tuple.h"
#include <string>

namespace mu2e {

struct CRVPhotonYieldConfig {
  using Name = fhicl::Name;
  using Comment = fhicl::Comment;

  fhicl::Atom<int> verbose{Name("verbose"), Comment("verbosity: 0 or 1")};
  fhicl::Atom<bool> useDb{Name("useDb"), Comment("use database or fcl")};
};

}  // namespace mu2e

#endif
