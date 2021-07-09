#include "Offline/MCDataProducts/inc/StageParticle.hh"

namespace mu2e {

  std::ostream& operator<<(std::ostream& os, const StageParticle& s) {
    return os<<"StageParticle("
             <<"parentId="<<s.parent()->id().asInt()
             <<", creationCode="<<s.creationCode().name()
             <<", pdgId="<<s.pdgId()
             <<", time="<<s.time()
             <<", position="<<s.position()
             <<", momentum="<<s.momentum()
             <<")";
  }

  std::ostream& operator<<(std::ostream& os, const StageParticleCollection& c) {
    for(const auto& s: c) {
      os<<s<<", ";
    }
    return os;
  }

}
