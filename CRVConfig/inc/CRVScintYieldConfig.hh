#ifndef STMConfig_CRVScintYieldConfig_hh
#define STMConfig_CRVScintYieldConfig_hh
//
// Fcl stanza for CRV ScintYield (scintillation yield spread of CRV counters) prodition
//
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Tuple.h"
#include <string>

namespace mu2e {

struct CRVScintYieldConfig {
  using Name = fhicl::Name;
  using Comment = fhicl::Comment;

  fhicl::Atom<int> verbose{Name("verbose"), Comment("verbosity: 0 or 1")};
  fhicl::Atom<bool> useDb{Name("useDb"), Comment("use database or fcl")};
};

}  // namespace mu2e

#endif
