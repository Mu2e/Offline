// Details for KalSeedMC-relevant structures
// Ed Callaghan, 2024

#include <Offline/MCDataProducts/inc/KalSeedMC.hh>

namespace mu2e{
  std::string const& TrkStrawHitProvenanceDetail::typeName(){
    static const std::string rv = "Provenance";
    return rv;
  };

  static const std::map<TrkStrawHitProvenanceDetail::enum_type, std::string> nam{
    std::make_pair(TrkStrawHitProvenance::unknown,    "unknown"),
    std::make_pair(TrkStrawHitProvenance::Simulation, "Simulation"),
    std::make_pair(TrkStrawHitProvenance::Mixed,      "Mixed"),
    std::make_pair(TrkStrawHitProvenance::External,   "External")
  };

  std::map<TrkStrawHitProvenance::enum_type, std::string> const& TrkStrawHitProvenanceDetail::names(){
    return nam;
  };
}
