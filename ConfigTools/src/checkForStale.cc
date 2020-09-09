//
// Given a list of keys and a SimpleConfig object, count how
// many of the keys have a value of true. Throw if more than
// one is true.  Optionally, throw if none are true.
//
//
// Original author Rob Kutschke
//

// Framework includes
#include "canvas/Utilities/Exception.h"

// Mu2e includes
#include "ConfigTools/inc/checkForStale.hh"
#include "ConfigTools/inc/SimpleConfig.hh"

using namespace std;

namespace mu2e {

  void checkForStale ( const vector<string>& staleNames,
                       const SimpleConfig&             config,
                       bool  throwOnMatch ){
    
    // Count how many of the matches there are

    int count(0);
    vector<string> variables;
    config.getNames( variables );

    for_each( variables.begin(),
              variables.end(),
              [&](string var) {
                if ( find( staleNames.begin(), staleNames.end(), var) != variables.end() ) {
                  cout << " Deprecated variable used: " << var << endl;
                  count++;
                }
              } );

    if ( count > 0 && throwOnMatch ) 
      throw cet::exception("GEOM") <<
        " Deprecated variables used in geom_01.txt " ;

  }

  void checkForStale ( string staleString,
                       const SimpleConfig& config,
                       bool  throwOnMatch ){
    
    // Count how many of the matches there are

    int count(0);
    vector<string> variables;
    config.getNames( variables );

    for_each( variables.begin(),
              variables.end(),
              [&](string var) {
                if ( var.find( staleString ) != string::npos ) {
                  cout << " Deprecated variable used: " << var << endl;
                  count++;
                }
              } );

    if ( count > 0 && throwOnMatch ) 
      throw cet::exception("GEOM") <<
        " Deprecated variables used in geom_01.txt " ;

  }

} // end namespace mu2e
