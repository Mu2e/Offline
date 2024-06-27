// Flag to designate a digi as having been produced from simulation,
// read from preexisting data, or a hybridization fo the two
// Ed Callaghan, 2024

#include <Offline/MCDataProducts/inc/DigiProvenance.hh>

namespace mu2e{
  std::string const& DigiProvenanceDetail::typeName(){
    static const std::string rv = "Provenance";
    return rv;
  }

  static const std::map<DigiProvenanceDetail::enum_type, std::string> nam{
    std::make_pair(DigiProvenance::unknown,    "unknown"),
    std::make_pair(DigiProvenance::Simulation, "Simulation"),
    std::make_pair(DigiProvenance::Mixed,      "Mixed"),
    std::make_pair(DigiProvenance::External,   "External")
  };

  std::map<DigiProvenance::enum_type, std::string> const& DigiProvenanceDetail::names(){
    return nam;
  }

  DigiProvenance::DigiProvenance(StringedDigiProvenance provenance)
      : StringedDigiProvenance(provenance){
        /**/
    }

  bool DigiProvenance::ContainsSimulation(){
    DigiProvenanceDetail::enum_type id = this->id();
    bool rv = ((id == DigiProvenanceDetail::Simulation) || (id == DigiProvenanceDetail::Simulation));
    return rv;
  }
}
