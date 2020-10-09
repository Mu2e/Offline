//
// Enum-matched-to-String class to identify the species
// of a reconstructable track.  Here a species is both
// a particle type (e/mu/pi/K/p) and a charge.
//
// See header file for more details.
//
// Contact person, Rob Kutschke
//

#include <type_traits>
#include <utility>

#include "Mu2eUtilities/inc/TrkSpecies.hh"

namespace mu2e {

  std::string const& TrkSpeciesDetail::typeName() {
    static const std::string type("TrkSpecies");
    return type;
  }

  std::map<TrkSpecies::enum_type,std::string> const& TrkSpeciesDetail::names(){
    static const std::map<TrkSpecies::enum_type,std::string> nam{
      std::make_pair(TrkSpecies::unknown,      "unknown"   ),
      std::make_pair(TrkSpecies::e_minus,      "eMinus"    ),
      std::make_pair(TrkSpecies::e_plus,       "ePlus"     ),
      std::make_pair(TrkSpecies::mu_minus,     "muMinus"   ),
      std::make_pair(TrkSpecies::mu_plus,      "muPlus"    ),
      std::make_pair(TrkSpecies::pi_plus,      "piPlus"    ),
      std::make_pair(TrkSpecies::pi_minus,     "piMinus"   ),
      std::make_pair(TrkSpecies::K_plus,       "KPlus"     ),
      std::make_pair(TrkSpecies::K_minus,      "KMinus"    ),
      std::make_pair(TrkSpecies::p_plus,       "pPlus"     ),
      std::make_pair(TrkSpecies::anti_p_minus, "pMinus"    ),
      std::make_pair(TrkSpecies::deuterium,    "deuterium" ),
      std::make_pair(TrkSpecies::tritium,      "tritium"   ),
      std::make_pair(TrkSpecies::He3,          "He3"       ),
      std::make_pair(TrkSpecies::He4,          "He4"       )
    };
    return nam;
  }

}
