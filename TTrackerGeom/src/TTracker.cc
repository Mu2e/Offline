//
// Geometry and identifier info about an TTracker.
//
//
// $Id: TTracker.cc,v 1.7 2011/05/26 22:08:09 genser Exp $
// $Author: genser $
// $Date: 2011/05/26 22:08:09 $
//
// Original author Rob Kutschke
//

#include "TTrackerGeom/inc/TTracker.hh"

using namespace std;

namespace mu2e {

  void TTracker::fillPointers () const{
    for ( size_t i=0; i<_devices.size(); ++i){
      _devices[i].fillPointers(*this);
    }
  }

  // Envelope that holds one device ("TTrackerDeviceEnvelope")
  TubsParams TTracker::getDeviceEnvelopeParams() const{

    double halfThick = _supportParams.halfThickness() + 2.*_manifoldHalfLengths[2];
    return TubsParams( _envelopeInnerRadius,
                       _supportParams.outerRadius(),
                       halfThick);
  }

  // Envelope that holds the full TTracker ("TrackerMother")
  TubsParams TTracker::getTrackerEnvelopeParams() const{

    // Envelope of a single device.
    TubsParams deviceEnvelope = getDeviceEnvelopeParams();

    // Full length from center to center of the first and last devices.
    double fullLength = _devices.back().origin().z()-_devices.front().origin().z();

    // Remember the thickness of the devices.
    double halfLength = fullLength/2. + deviceEnvelope.zHalfLength();

    return TubsParams( deviceEnvelope.innerRadius(),
                       deviceEnvelope.outerRadius(),
                       halfLength);
  }


} // namespace mu2e
