//
// Given a list of keys and a SimpleConfig object, count how
// many of the keys have a value of true. Throw if more than
// one is true.  Optionally, throw if none are true.
//
// $Id: requireUniqueKey.cc,v 1.1 2012/07/15 22:00:36 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/07/15 22:00:36 $
//
// Original author Rob Kutschke
//

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Utilities/Exception.h"

// Mu2e includes
#include "ConfigTools/inc/requireUniqueKey.hh"
#include "ConfigTools/inc/SimpleConfig.hh"

using namespace std;

namespace mu2e {

  int requireUniqueKey ( const std::vector<std::string>& keys,
                         const SimpleConfig&             config,
                         bool  throwOnZero ){

    // Count how many of the keys are true?
    int count(0);
    string found;
    for ( size_t i=0; i<keys.size(); ++i ){
      if ( config.getBool(keys[i],false) ) {
        if ( !found.empty() ) found += " ";
        found += keys[i];
        ++count;
      }
    }

    // Throw if an error detected.
    if ( count == 0 ) {

      if ( throwOnZero ){
        cet::exception exception("CONFIG");
        exception << "None of the configuration parameters is true: ";
        for ( size_t i=0; i<keys.size(); ++i){
          if ( i != 0 ) exception << " ";
          exception << keys[i];
        }
        exception << "\n";
        throw exception;
      }

    } else if ( count > 1 ) {
        throw cet::exception("CONFIG")
          << "More than one of the requested configuration parameters is true: "
          << found
          << "\n";
    }

    // Successful return.
    return count;
  }

} // end namespace mu2e
