//
// This is a place to put additional information produced by HitMaker,
//
// $Id: ExtMonUCITofHitMCTruth.cc,v 1.1 2011/12/30 20:31:46 youzy Exp $
// $Author: youzy $
// $Date: 2011/12/30 20:31:46 $
//
//

// C++ includes
#include <ostream>

// Framework includes.
#include "cetlib_except/exception.h"

// Mu2e includes
#include "MCDataProducts/inc/ExtMonUCITofHitMCTruth.hh"

using namespace std;

namespace mu2e {

  // Print the information found in this hit.
  void ExtMonUCITofHitMCTruth::print( ostream& ost, bool doEndl ) const {

    ost << "ExtMonTof Hit MC:"
        << " Station id: " << _stationId
        << " Segment id: " << _segmentId
        << " time: "       << _time
        << " energyDep: "  << _energyDep
        << " track id: "   << _trackId
        << " pdg id: "     << _pdgId
        << " position "    << _position
        << " momentum "    << _momentum
        << " vertex "      << _vertex
        << " vertexMomentum " << _vertexMomentum
        << " vertexTime "     << _vertexTime
        << " isPrimary "      << _isPrimary
        << " org track id:  " << _orgTrackId
        << " org pdg id: "    << _orgPdgId
        << " org vertex "     << _orgVertex
        << " org vertex momentum " << _orgVertexMomentum
        << " org time " << _orgTime;

    if ( doEndl ){
      ost << endl;
    }

  }

} // namespace mu2e
