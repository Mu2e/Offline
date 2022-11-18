#ifndef STMConfig_CRVOrdinalConfig_hh
#define STMConfig_CRVOrdinalConfig_hh
//
// Fcl stanza for CRV Ordinal (ROC numbering) prodition
//
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Tuple.h"
#include <string>

namespace mu2e {

struct CRVOrdinalConfig {
  using Name = fhicl::Name;
  using Comment = fhicl::Comment;

  fhicl::Atom<int> verbose{Name("verbose"), Comment("verbosity: 0 or 1")};
  fhicl::Atom<bool> useDb{Name("useDb"), Comment("use database or fcl")};

};

}  // namespace mu2e

#endif
