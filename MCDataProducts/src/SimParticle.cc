#include "MCDataProducts/inc/SimParticle.hh"

namespace mu2e {
  std::vector<SimParticle::key_type>  SimParticle::daughterIds() const {
    std::vector<key_type> res;
    res.reserve(_daughterSims.size());
    for(const auto i: _daughterSims) {
      res.push_back( key_type(i.key()) );
    }
    return res;
  }
}
