// Flag to designate a digi as having been produced from simulation,
// read from preexisting data, or a hybridization fo the two
// Ed Callaghan, 2024

#ifndef MCDataProducts_DigiProvenance_hh
#define MCDataProducts_DigiProvenance_hh

// stl
#include <map>
#include <string>

// mu2e
#include "Offline/GeneralUtilities/inc/EnumToStringSparse.hh"

namespace mu2e{
  // enum equipped with std::string descriptions
  class DigiProvenanceDetail{
    public:
      enum enum_type {unknown=0, Simulation, Mixed, External};
      static std::string const& typeName();
      static std::map<enum_type, std::string> const& names();
  };
  using DigiProvenance = EnumToStringSparse<DigiProvenanceDetail>;

  bool containsSimulation(const DigiProvenance&);
}

#endif
