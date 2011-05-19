//
// Geometry and identifier info about an TTracker.
//
//
// $Id: TTracker.cc,v 1.5 2011/05/19 22:23:06 wb Exp $
// $Author: wb $
// $Date: 2011/05/19 22:23:06 $
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

    // Padding is half the size of the padding for the envelope of the full tracker.
    static double const pad = 0.00005;
    double halfThick = _supportParams.halfThickness() + 2.*_manifoldHalfLengths[2];
    return TubsParams( _envelopeInnerRadius-pad,
                       _supportParams.outerRadius()+pad,
                       halfThick+pad);
  }

  // Envelope that holds the full TTracker ("TrackerMother")
  TubsParams TTracker::getTrackerEnvelopeParams() const{

    // Envelope of a single device.
    TubsParams deviceEnvelope = getDeviceEnvelopeParams();

    // Full length from center to center of the first and last devices.
    double fullLength = _devices.back().origin().z()-_devices.front().origin().z();

    // Remember the thickness of the devices.
    double halfLength = fullLength/2. + deviceEnvelope.zHalfLength;

    static double const pad = 0.0001;
    return TubsParams( deviceEnvelope.innerRadius-pad,
                       deviceEnvelope.outerRadius+pad,
                       halfLength+pad);

    ;
  }


} // namespace mu2e
