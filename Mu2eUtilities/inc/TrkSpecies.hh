#ifndef Mu2eUtilities_TrkSpecies_hh
#define Mu2eUtilities_TrkSpecies_hh
//
// Enum-matched-to-String class to identify the species
// of a reconstructable track.  Here a species is both
// a particle type (e/mu/pi/K/p) and a charge.
//
// This is the same information provided by the class
//    $BTRK_INC/BTrk/TrkBase/TrkParticle.hh
// Someday we should integrate the two.
//
// I could have called this mu2e::TrkParticle to emphasize
// the connection to ::TrkParticle but that seemd too dangerous:
// because the BTRK name is in the global namespace, there is some
// code might compile against the wrong include file.
//
// Contact person, Rob Kutschke
//

#include <map>
#include <string>

#include "GeneralUtilities/inc/EnumToStringSparse.hh"

namespace mu2e {

  class TrkSpeciesDetail{
  public:

    enum enum_type {
      unknown      =     0,
      e_minus      =    11,
      e_plus       =   -11,
      mu_minus     =    13,
      mu_plus      =   -13,
      pi_plus      =   211,
      pi_minus     =  -211,
      K_plus       =   321,
      K_minus      =  -321,
      p_plus       =  2212,
      anti_p_minus = -2212,
      deuterium    =  1000010020,
      tritium      =  1000010030,
      He3          =  1000020030,
      He4          =  1000020040
    };

    static std::string const& typeName();

    static std::map<enum_type,std::string> const& names();

  };

  typedef EnumToStringSparse<TrkSpeciesDetail> TrkSpecies;
}

#endif
