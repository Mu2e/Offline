// Geometrical info for the extinction monitor
//
// Andrei Gaponenko, 2011

#ifndef EXTMONFNAL_HH
#define EXTMONFNAL_HH

#include <vector>

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/RotationX.h"

#include "canvas/Persistency/Common/Wrapper.h"

#include "Offline/Mu2eInterfaces/inc/Detector.hh"

#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALPixelChip.hh"
#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALModule.hh"
#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALPlane.hh"
#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALPlaneStack.hh"
#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALMagnet.hh"

namespace mu2e {

  class ExtMonFNALPixelId;

  namespace ExtMonFNAL {

    class ExtMonMaker;

    class ExtMon : virtual public Detector {

    public:

      // all modules are the same
      const ExtMonFNALModule& module() const { return module_; }

      // all chips are the same
      const ExtMonFNALPixelChip& chip() const { return chip_; }

      const ExtMonFNALPlaneStack& up() const { return up_; }
      const ExtMonFNALPlaneStack& dn() const { return dn_; }

      const ExtMonFNALPlane& plane(int p) const {
        if (p < dn_.size())
          return dn_.planes()[p];
        else return up_.planes()[p - dn_.size()];
      }


      const ExtMonFNALMagnet& spectrometerMagnet() const { return spectrometerMagnet_; }

      // Location of the detector == that of the upstream stack.
      CLHEP::Hep3Vector detectorCenterInMu2e() const;
      const CLHEP::HepRotation& detectorRotationInMu2e() const;

      const std::vector<double>& detectorMotherHS() const { return detectorMotherHS_; }
      const CLHEP::Hep3Vector detectorMotherCenterInMu2e() const { return detectorMotherCenterInMu2e_; }
      const CLHEP::HepRotation& detectorMotherRotationInMu2e() const { return spectrometerMagnet().magnetRotationInMu2e(); }

      // Coordinate conversion to/from the Mu2e frame
      // The ExtMonFNAL frame is defined in the following way:
      //
      // - The (0,0,0) is the reference point of the upstream plane stack,
      //
      // - The z_em axis is along the upstream stack axis (perpendicular to the plane planes)
      // - The x_em axis in in the horizontal plane
      // - The y_em axis forms a right-handed (x_em, y_em, z_em) frame
      //
      // The rotation of the ExtMon system w.r.t. Mu2e is "small",
      // that is, the projection of z_em on z_mu2e is positive and
      // similar for x, y.

      CLHEP::Hep3Vector mu2eToExtMon_position(const CLHEP::Hep3Vector& mu2epos) const;
      CLHEP::Hep3Vector mu2eToExtMon_momentum(const CLHEP::Hep3Vector& mu2emom) const;

      CLHEP::Hep3Vector extMonToMu2e_position(const CLHEP::Hep3Vector& pos) const;
      CLHEP::Hep3Vector extMonToMu2e_momentum(const CLHEP::Hep3Vector& mom) const;

      //----------------------------------------------------------------
      CLHEP::Hep3Vector upStackToExtMon_position(const CLHEP::Hep3Vector& uppos) const {
        return uppos;
      }
      CLHEP::Hep3Vector dnStackToExtMon_position(const CLHEP::Hep3Vector& dnpos) const {
        return dnToExtMonCoordinateRotation_ * dnpos;
      }
      CLHEP::Hep3Vector stackToExtMon_position(const CLHEP::Hep3Vector& pos) const {
        return (pos.z() < 0) ? dnStackToExtMon_position(pos) : upStackToExtMon_position(pos);
      }

      //----------------------------------------------------------------
      // Pixel center in the coordinate system of its PlaneStack

      CLHEP::Hep3Vector pixelPositionInPlaneStack(const ExtMonFNALPixelId& id) const;

      //----------------------------------------------------------------
      // Redundant convenience accessors
      unsigned int nplanes() const { return dn_.nplanes() + up_.nplanes(); }

      unsigned int nmodules() const {
        unsigned int nmod = 0;
        for (unsigned iplane = 0; iplane < up_.nplanes(); iplane++)
          nmod += up_.planes()[iplane].nModules();
        for (unsigned iplane = 0; iplane < dn_.nplanes(); iplane++)
          nmod += dn_.planes()[iplane].nModules();
        return nmod;
      }
      CLHEP::Hep3Vector planeCenterInExtMon(unsigned iplane) const;

      bool samePlaneStack(unsigned plane1, unsigned plane2) {
        bool dn1 = (plane1 < dn_.nplanes());
        bool dn2 = (plane2 < dn_.nplanes());
        return !(dn1^dn2);
      }

    private:
      friend class ExtMonMaker;
      ExtMon() {}

      // For persistency
      template<class T> friend class art::Wrapper;

      ExtMonFNALPixelChip chip_;
      ExtMonFNALModule module_;
      ExtMonFNALPlaneStack up_;
      ExtMonFNALPlaneStack dn_;
      ExtMonFNALMagnet spectrometerMagnet_;
      CLHEP::HepRotationX dnToExtMonCoordinateRotation_;
      std::vector<double> detectorMotherHS_;
      CLHEP::Hep3Vector detectorMotherCenterInMu2e_;
    };

    //================================================================

  }
}

#endif/*EXTMONFNAL_HH*/
