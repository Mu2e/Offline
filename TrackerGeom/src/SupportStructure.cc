//
// A model for the Tracker supports described in Mu2e-doc-888-v5 (current as of Dec, 2012).
//
//  $Id: SupportStructure.cc,v 1.3 2014/01/06 20:50:40 kutschke Exp $
//  $Author: kutschke $
//  $Date: 2014/01/06 20:50:40 $
//
//  Original author Rob Kutschke
//

#include "TrackerGeom/inc/SupportStructure.hh"
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
    _stiffRings(),
    _beamBody(),
    _beamServices(){
  }

  void SupportStructure::print ( std::ostream& ost ) const{

    ost << "\nNumber of Stiffening rings: " << _stiffRings.size() << std::endl;
    for ( size_t i=0; i<_stiffRings.size(); ++i ){
      ost << "\nStiff ring: "
          << i << "    "
          << _stiffRings[i]
          << std::endl;
    }


    for ( size_t i=0; i<_beamBody.size(); ++i ){
      ost << "\nBeam body: "
          << i << "    "
          << _beamBody[i]
          << std::endl;
    }

    for ( size_t i=0; i<_beamServices.size(); ++i ){
      ost << "\nBeam services: "
          << i << "    "
          << _beamServices[i]
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
