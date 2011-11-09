// Geometrical info for the extinction monitor
// 
// Andrei Gaponenko, 2011

#ifndef EXTMONFNAL_HH
#define EXTMONFNAL_HH

#include <vector>

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"

#include "GeometryService/inc/Detector.hh"

namespace mu2e {
  namespace ExtMonFNAL {
    
    class ExtMonMaker;

    class ExtMon : public Detector {
    public: 
      // implement Detector's method
      virtual std::string name() const { return "ExtMonFNAL"; }
      
      // 
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

      // Coordinate conversion to/from the Mu2e frame
      // The ExtMonFNAL frame is defined in the following way:
      // 
      // - The z_em axis is along the channel (i.e. the secondaries we 
      //   observe have large pz and small px, py).
      // - The x_em axis in in the horizontal plane
      // - The y_em axis forms a right-handed (x_em, y_em, z_em) frame.
      // 
      // - The (0,0,0) point is in the middle of the ExtMon detector.
      // 
      // In G4 this is the coordinate system of the ExtMonFNAL volume.
      CLHEP::Hep3Vector extMonToMu2ePoint( CLHEP::Hep3Vector const& v ) const;      
      CLHEP::Hep3Vector extMonToMu2eMomentum( CLHEP::Hep3Vector const& v ) const;      
      CLHEP::Hep3Vector mu2eToExtMonPoint( CLHEP::Hep3Vector const& v ) const;
      CLHEP::Hep3Vector mu2eToExtMonMomentum( CLHEP::Hep3Vector const& v ) const;

      // G4 implementation details
      const std::vector<double>&  logicalEnclosureHalfDim() const { return m_logicalEnclosureHalfDim; }
      CLHEP::Hep3Vector  const& offsetInParent() const { return m_offsetInParent; }
      CLHEP::HepRotation const& rotationInParent() const { return m_rotationInParent; }

      //----------------------------------------------------------------
    private: 
      friend class ExtMonMaker;
      // Private ctr: the class should be only obtained via ExtMonFNAL::ExtMonMaker.
      ExtMon(const std::vector<double>& logicalEnclosureHalfDim, 
	     const CLHEP::Hep3Vector& offsetInParent);

      std::vector<double> m_logicalEnclosureHalfDim;
      CLHEP::Hep3Vector m_offsetInParent;
      CLHEP::HepRotation m_rotationInParent;

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
    };

  }
}

#endif/*EXTMONFNAL_HH*/
