#include "ServicesGeom/inc/Pipe.hh"
#include <vector>

namespace mu2e {

  // genreflex persistency requires default ctr
  Pipe::Pipe() {}


  std::ostream& operator<<(std::ostream& os, const Pipe& ens) {
    os<<"Pipe("
      <<"Contains " << ens.getNPipes().size()
      <<" Types of pipes "
      <<" )";
    return os;
  }

} // namespace mu2e
