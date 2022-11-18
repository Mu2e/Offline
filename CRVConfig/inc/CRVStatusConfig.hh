#ifndef STMConfig_CRVStatusConfig_hh
#define STMConfig_CRVStatusConfig_hh
//
// Fcl stanza for CRV Status (bad channels) prodition
//
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Tuple.h"
#include <string>

namespace mu2e {

struct CRVStatusConfig {
  using Name = fhicl::Name;
  using Comment = fhicl::Comment;

  fhicl::Atom<int> verbose{Name("verbose"), Comment("verbosity: 0 or 1")};
  fhicl::Atom<bool> useDb{Name("useDb"), Comment("use database or fcl")};
  fhicl::Sequence<fhicl::Tuple<std::size_t, int>> status{
      Name("status"), Comment("sequence of tuples of channel and status_flag")};
};

} // namespace mu2e

#endif
