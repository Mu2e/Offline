// Geometrical info for the extinction monitor
//
// Andrei Gaponenko, 2011

#ifndef EXTMONFNAL_HH
#define EXTMONFNAL_HH

#include <vector>

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"

#include "art/Persistency/Common/Wrapper.h"

#include "Mu2eInterfaces/inc/Detector.hh"

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALSensorStack.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALMagnet.hh"

namespace mu2e {
  namespace ExtMonFNAL {

    class ExtMonMaker;

    class ExtMon : virtual public Detector {
    public:

      const ExtMonFNALSensorStack& up() const { return up_; }
      const ExtMonFNALSensorStack& dn() const { return dn_; }

      const ExtMonFNALMagnet& spectrometerMagnet() const { return spectrometerMagnet_; }

      // Location of the detector == that of the upstream stack.
      CLHEP::Hep3Vector detectorCenterInMu2e() const;
      const CLHEP::HepRotation& detectorRotationInMu2e() const;

      // Coordinate conversion to/from the Mu2e frame
      // The ExtMonFNAL frame is defined in the following way:
      //
      // - The (0,0,0) is the reference point of the upstream sensor stack,
      //
      // - The z_em axis is along the upstream stack axis (perpendicular to the sensor planes)
      // - The x_em axis in in the horizontal plane
      // - The y_em axis forms a right-handed (x_em, y_em, z_em) frame
      //
      // The rotation of the ExtMon system w.r.t. Mu2e is "small",
      // that is, the projection of z_em on z_mu2e is positive and
      // similar for x, y.
      //
      // In G4 this is the coordinate system of the ExtMonFNAL volume.

      CLHEP::Hep3Vector mu2eToExtMon_position(const CLHEP::Hep3Vector& mu2epos) const;
      CLHEP::Hep3Vector mu2eToExtMon_momentum(const CLHEP::Hep3Vector& mu2emom) const;

      //----------------------------------------------------------------
    private:
      friend class ExtMonMaker;
      ExtMon() {}

      // For persistency
      template<class T> friend class art::Wrapper;

      ExtMonFNALSensorStack up_;
      ExtMonFNALSensorStack dn_;
      ExtMonFNALMagnet spectrometerMagnet_;
    };

    //================================================================

  }
}

#endif/*EXTMONFNAL_HH*/
