#include "Offline/Mu2eKinKal/inc/StrawHitUpdaters.hh"
namespace mu2e {
  std::vector<std::string> StrawHitUpdaters::names_{"None","CADSHU","DriftANNSHU","BkgANNSHU","","","","","","PanelDiagSHU","Chi2SHU"};
  std::string const& StrawHitUpdaters::name(algorithm alg) {
    // 'unknown' is -1, so it must never be used as an index into names_
    static const std::string unknownname("unknown");
    auto ialg = static_cast<size_t>(alg);
    if(alg == unknown || ialg >= names_.size()) return unknownname;
    return names_[ialg];
  }
  StrawHitUpdaters::algorithm StrawHitUpdaters::algo(std::string const& name) {
    for(size_t alg=0; alg < StrawHitUpdaters::nalgos;++alg){
      if(names_[alg].compare(name)==0)
        return static_cast<StrawHitUpdaters::algorithm>(alg);
    }
    return StrawHitUpdaters::unknown;
  }
}
