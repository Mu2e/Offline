//
// Master function for checking for a valid geometry.
//
#include "Offline/Mu2eG4/inc/validGeometryOrThrow.hh"
#include "Offline/Mu2eG4/inc/validPolyCones.hh"

#include "cetlib_except/exception.h"


void mu2e::validGeometryOrThrow( int printLevel ){

  // Calls to individual tests; add more here as needed.
  bool pConesOK = validPolyCones( printLevel );

  // Throw an exception if any test fails.
  // Put detailed printout in the indvidual test functions.
  if ( !pConesOK ){
    throw cet::exception("GEOM")
      << "Mu2eG4::validGeometryOrThrow Error: inconsistent geometry\n";
  }

}
