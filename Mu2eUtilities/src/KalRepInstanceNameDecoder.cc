//
// Decode the instance name of a KalRepCollection
// or a KalRepPtrCollection.
//
// Contact person, Rob Kutschke
//

#include "GeneralUtilities/inc/EnumToStringSparse.hh"
#include "cetlib_except/exception.h"

#include "Mu2eUtilities/inc/KalRepInstanceNameDecoder.hh"

mu2e::KalRepInstanceNameDecoder::KalRepInstanceNameDecoder( std::string const& instanceName ):
  instanceName_(instanceName),
  direction_(),
  particleType_(){

  // Particle species field of the instance name; includes the particle name and charge.
  std::string species;

  // Identify direction; extract species field.
  if ( instanceName_.substr(0,8) == "Upstream" ){
    direction_ = TrkFitDirection::upstream;
    species = instanceName_.substr(8);
  } else if ( instanceName_.substr(0,10) == "Downstream" ){
    direction_ = TrkFitDirection::downstream;
    species = instanceName_.substr(10);
  } else{
    throw cet::exception("TRACKERINSTANCE")
      <<"Cannot parse direction from track data product instance name = " << instanceName <<"\n";
  }

  particleType_ = TrkSpecies(species);
  if ( particleType_ == TrkSpecies::unknown ){
    throw cet::exception("TRACKERINSTANCE")
      <<"Cannot parse species from track data product instance name = " << instanceName <<"\n";
  }
  if ( species.find("Plus") != std::string::npos ) {
    charge_ = +1;
  } else if ( species.find("Minus") != std::string::npos ){
    charge_ = -1;
  } else{
    throw cet::exception("TRACKERINSTANCE")
      <<"Cannot parse charge from track data product instance name = " << instanceName <<"\n";
  }

} // end c'tor mu2e::KalRepInstanceNameDecoder
