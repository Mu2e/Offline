//
// Enum-matched-to-String class for defining which version of the Tracker
// support structure is to be used.
//
//   $Id: SupportModel.cc,v 1.1 2013/01/07 03:54:28 kutschke Exp $
//   $Author: kutschke $
//   $Date: 2013/01/07 03:54:28 $
//
// Contact person, Rob Kutschke
//

#include "TrackerGeom/inc/SupportModel.hh"

namespace mu2e {

  std::string const& SupportModelDetail::typeName() {
    static std::string type("SupportModel");
    return type;
  }

  std::map<SupportModelDetail::enum_type,std::string> const& SupportModelDetail::names(){

    static std::map<enum_type,std::string> nam;

    if ( nam.empty() ){
      nam[unknown]    = "unknown";
      nam[simple]     = "simple";
      nam[detailedv0] = "detailedv0";
    }

    return nam;
  }

}
