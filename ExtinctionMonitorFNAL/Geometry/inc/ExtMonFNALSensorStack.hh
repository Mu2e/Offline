// Stack of the silicon sensor planes.
//
// Andrei Gaponenko, 2012

#ifndef EXTMONFNALSENSORSTACK_HH
#define EXTMONFNALSENSORSTACK_HH

#include <vector>

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"

#include "art/Persistency/Common/Wrapper.h"

namespace mu2e {

  namespace ExtMonFNAL { class ExtMonMaker; }
  namespace ExtMonFNAL { class ExtMon; }

  class ExtMonFNALSensorStack {
  public:

    unsigned nplanes() const { return m_sensor_zoffset.size(); }
    const std::vector<double>& sensor_zoffset() const { return m_sensor_zoffset; }
    const std::vector<double>& sensor_xoffset() const { return m_sensor_xoffset; }
    const std::vector<double>& sensor_yoffset() const { return m_sensor_yoffset; }

    const std::vector<double>& readout_halfdz() const { return m_readout_halfdz; }

    // offset of sensor center wrt the ref point
    CLHEP::Hep3Vector   sensorOffsetInStack(unsigned iplane) const;

    // Positioning of the stack: reference point and rotation.
    CLHEP::Hep3Vector refPointInMu2e() const { return m_stackRefPointInMu2e; }
    CLHEP::HepRotation const& rotationInMu2e() const { return m_stackRotationInMu2e; }

    // Test materials for MARS activation studies
    const std::vector<double>& testMaterialHalfSize() const { return m_testMaterialHalfSize; }
    double distanceToTestMaterials() const { return distanceToTestMaterials_; }
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

    CLHEP::Hep3Vector mu2eToStack_position(const CLHEP::Hep3Vector& mu2epos) const;
    CLHEP::Hep3Vector mu2eToStack_momentum(const CLHEP::Hep3Vector& mu2emom) const;

    //----------------------------------------------------------------
  private:
    ExtMonFNALSensorStack();
    friend class ExtMonFNAL::ExtMon;
    friend class ExtMonFNAL::ExtMonMaker;

    // For persistency
    template<class T> friend class art::Wrapper;

    CLHEP::HepRotation m_stackRotationInMu2e;
    CLHEP::Hep3Vector m_stackRefPointInMu2e;

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
    double distanceToTestMaterials_;
    double m_testMaterialPitch;
    std::vector<std::string> m_testMaterialNames;

    // data for coordinate system transformations: inverse of stack rotation
    CLHEP::HepRotation m_coordinateRotationInMu2e;
  };

  //================================================================

} // namespace mu2e

#endif/*EXTMONFNALSENSORSTACK_HH*/
