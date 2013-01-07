//
// A model for the TTracker supports described in Mu2e-doc-888-v5 (current as of Dec, 2012).
//
//  $Id: SupportStructure.cc,v 1.1 2013/01/07 03:56:08 kutschke Exp $
//  $Author: kutschke $
//  $Date: 2013/01/07 03:56:08 $
//
//  Original author Rob Kutschke
//

#include "TTrackerGeom/inc/SupportStructure.hh"
#include <iostream>

namespace mu2e {

  SupportStructure::SupportStructure():
    _centerPlate(),
    _gasUpstream(),
    _gasDownstream(),
    _innerRing(),
    _outerRingUpstream(),
    _outerRingDownstream(),
    _coverUpstream(),
    _coverDownstream(),
    _innerChannelUpstream(),
    _innerChannelDownstream(),
    _g10Upstream(),
    _g10Downstream(),
    _cuUpstream(),
    _cuDownstream(),
    _endRingUpstream(),
    _endRingDownstream(),
    _staveBody(),
    _staveServices(){
  }

  void SupportStructure::print ( std::ostream& ost ) const{

    ost << "\nEnd Ring Upstream:      " <<  _endRingUpstream       << std::endl;
    ost << "\nEnd Ring Downstream:    " <<  _endRingDownstream     << std::endl;

    for ( size_t i=0; i<_staveBody.size(); ++i ){
      ost << "\mStave body: "
          << i << "    "
          << _staveBody[i]
          << std::endl;
    }

    for ( size_t i=0; i<_staveServices.size(); ++i ){
      ost << "\nStave services: "
          << i << "    "
          << _staveServices[i]
          << std::endl;
    }

    ost << "\nCenterPlate:            " << _centerPlate            << std::endl;
    ost << "\nGas Upstream:           " << _gasUpstream            << std::endl;
    ost << "\nGas Downstream:         " << _gasDownstream          << std::endl;
    ost << "\nInner Ring:             " << _innerRing              << std::endl;
    ost << "\nOuterRing upstream:     " << _outerRingUpstream      << std::endl;
    ost << "\nOuterRign downstream:   " << _outerRingDownstream    << std::endl;
    ost << "\nCover upstream:         " << _coverUpstream          << std::endl;
    ost << "\nCover downstream:       " << _coverDownstream        << std::endl;
    ost << "\nChannel upstream:       " << _innerChannelUpstream   << std::endl;
    ost << "\nChannel downstream:     " << _innerChannelDownstream << std::endl;
    ost << "\nG10 upstream:           " <<  _g10Upstream           << std::endl;
    ost << "\nG10 downstream:         " <<  _g10Downstream         << std::endl;
    ost << "\nCu upstream:            " <<  _cuUpstream            << std::endl;
    ost << "\nCu downstream:          " <<  _cuDownstream          << std::endl;

  }

}
