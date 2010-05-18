//
// Geometry and identifier info about an TTracker.
//
//
// $Id: TTracker.cc,v 1.2 2010/05/18 21:16:51 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/05/18 21:16:51 $
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

  // Envelope that holds one device:
  TubsParams TTracker::getDeviceEnvelopeParams() const{

    // Padding is half the size of the padding for the envelope of the full tracker.
    double pad = 0.05;
    double halfThick = _supportParams.halfThickness + 2.*_manifoldHalfLengths[2];
    return TubsParams( _envelopeInnerRadius-pad,
                       _supportParams.outerRadius+pad,
                       halfThick+pad);
  }

  // Envelope that holds the full TTracker
  TubsParams TTracker::getTrackerEnvelopeParams() const{

    // Envelope of a single device.
    TubsParams deviceEnvelope = getDeviceEnvelopeParams();

    // Full length from center to center of the first and last devices.
    double fullLength = _devices.back().origin().z()-_devices.front().origin().z();

    // Remember the thickness of the devices.
    double halfLength = fullLength/2. + deviceEnvelope.zHalfLength;

    double pad = 0.1;
    return TubsParams( deviceEnvelope.innerRadius-pad,
                       deviceEnvelope.outerRadius+pad,
                       halfLength+pad);

    ;
  }
  

} // namespace mu2e
