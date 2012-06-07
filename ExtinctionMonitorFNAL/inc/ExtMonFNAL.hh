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

namespace mu2e {
  namespace ExtMonFNAL {

    class ExtMonMaker;

    class ExtMon : virtual public Detector {
    public:
      unsigned nplanes() const { return m_sensor_zoffset.size(); }
      const std::vector<double>& sensor_zoffset() const { return m_sensor_zoffset; }
      const std::vector<double>& sensor_xoffset() const { return m_sensor_xoffset; }
      const std::vector<double>& sensor_yoffset() const { return m_sensor_yoffset; }

      const std::vector<double>& sensor_halfdx() const { return m_sensor_halfdx; }
      const std::vector<double>& sensor_halfdy() const { return m_sensor_halfdy; }
      const std::vector<double>& sensor_halfdz() const { return m_sensor_halfdz; }

      const std::vector<double>& readout_halfdz() const { return m_readout_halfdz; }

      // same info as above cooked for nestBox()
      std::vector<double> sensorHalfSize(unsigned iplane) const;
      CLHEP::Hep3Vector sensorOffsetInParent(unsigned iplane) const;

      // Location of the detector
      CLHEP::Hep3Vector detectorCenterInMu2e() const { return m_detectorCenterInMu2e; }
      CLHEP::HepRotation const& detectorRotationInMu2e() const { return m_detectorRotationInMu2e; }
      // the size is computed from sensor pars above
      const std::vector<double>& detectorHalfSize() const { return m_detectorHalfSize; }

      // Test materials for MARS activation studies
      const std::vector<double>& testMaterialHalfSize() const { return m_testMaterialHalfSize; }
      double testMaterialDistanceToDetector() const { return m_testMaterialDistanceToDetector; }
      double testMaterialPitch() const { return m_testMaterialPitch; }
      const std::vector<std::string>&  testMaterialNames() const { return m_testMaterialNames; }

      // Coordinate conversion to/from the Mu2e frame
      // The ExtMonFNAL frame is defined in the following way:
      //
      // - The (0,0,0) point is in at the intersection of collimator2 axis with
      //   the surface of the ExtMonFNAL room wall.
      //
      // - The z_em axis is along the collimator2 channel
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
      // Private ctr: the class should be only obtained via ExtMonFNAL::ExtMonMaker.
      ExtMon() {}

      // Or read back from persistent storage
      template<class T> friend class art::Wrapper;

      std::vector<double> m_detectorHalfSize;
      CLHEP::HepRotation m_detectorRotationInMu2e;
      CLHEP::Hep3Vector m_detectorCenterInMu2e;

      // Sensor center positions
      std::vector<double> m_sensor_zoffset;
      std::vector<double> m_sensor_xoffset;
      std::vector<double> m_sensor_yoffset;

      // Sensor size
      std::vector<double> m_sensor_halfdx;
      std::vector<double> m_sensor_halfdy;
      std::vector<double> m_sensor_halfdz;

      // Readout electronics is (at the moment) created with the same (x,y) size as the sensor,
      // as a box parallel to the sensor.  One thing that remains to be specified is:
      std::vector<double> m_readout_halfdz;

      // Test materials for MARS activation studies
      std::vector<double> m_testMaterialHalfSize;
      double m_testMaterialDistanceToDetector;
      double m_testMaterialPitch;
      std::vector<std::string> m_testMaterialNames;

      // data for coordinate system transformations
      CLHEP::Hep3Vector m_coordinateCenterInMu2e;
      CLHEP::HepRotation m_coordinateRotationInMu2e;
    };

  }
}

#endif/*EXTMONFNAL_HH*/
