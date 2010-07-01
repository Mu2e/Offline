// 
// First version of a hit as described by Mu2e-doc-900.
//
// $Id: StrawHit.cc,v 1.1 2010/07/01 13:34:57 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/07/01 13:34:57 $
//
// Original author Rob Kutschke
//

// C++ includes
#include <ostream>

// Framework includes.
#include "FWCore/Utilities/interface/Exception.h"

// Mu2e includes
#include "ToyDP/inc/StrawHit.hh"

using namespace std;

namespace mu2e {

  // Print the information found in this hit.
  void StrawHit::print( ostream& ost, bool doEndl ) const {

    ost << "traw Hit:"
        << " idx: "      << _strawIndex
        << " time: "     << _time
        << " dt: "       << _dt
        << " eDep: "     << _energyDep
        << " pre: ";

    printPrecursorIndices(ost);

    if ( doEndl ){
      ost << endl;
    }
    
  }
  
  // Print the values of the elements of precursorIndices.
  void StrawHit::printPrecursorIndices( ostream& ost) const {

    ost << "(";
    for ( vector<uint32_t>::size_type i=0;
          i<_precursorIndices.size(); ++i ){
      if ( i > 0 ){
        ost << ",";
      }
      ost << _precursorIndices[i];
    }
    ost << ")";
  }

} // namespace mu2e
