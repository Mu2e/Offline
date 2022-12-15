#include "Offline/Mu2eKinKal/inc/StrawHitUpdaters.hh"
namespace mu2e {
  std::vector<std::string> StrawHitUpdaters::names_{"None","CADSHU","DriftANNSHU","BkgANNSHU","","","","","","","Chi2SHU"};
  std::string const& StrawHitUpdaters::name(algorithm alg) {
    return names_[static_cast<size_t>(alg)];
  }
  StrawHitUpdaters::algorithm StrawHitUpdaters::algo(std::string const& name) {
    for(size_t alg=0; alg < StrawHitUpdaters::nalgos;++alg){
      if(names_[alg].compare(name)==0)
        return static_cast<StrawHitUpdaters::algorithm>(alg);
    }
    return StrawHitUpdaters::unknown;
  }
}
