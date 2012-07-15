// Andrei Gaponenko, 2011

#include "ExtinctionMonitorFNAL/inc/ExtMonFNAL_Maker.hh"
#include "ExtinctionMonitorFNAL/inc/ExtMonFNAL.hh"
#include "ExtinctionMonitorFNAL/inc/ExtMonFNALBuilding.hh"

#include "cetlib/exception.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include "ConfigTools/inc/SimpleConfig.hh"

//#define AGDEBUG(stuff) std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<": "<<stuff<<std::endl;
#define AGDEBUG(stuff)

namespace mu2e {
  namespace ExtMonFNAL {

    std::auto_ptr<ExtMon> ExtMonMaker::make(const SimpleConfig& config, const ExtMonFNALBuilding& room)
    {
      AGDEBUG("ExtMonFNAL maker begin");
      const int verbose = config.getInt("extMonFNAL.verbosityLevel", 0);

      std::auto_ptr<ExtMon> det(new ExtMon());

      config.getVectorDouble("extMonFNAL.sensor_zoffset", det->m_sensor_zoffset, -1);
      config.getVectorDouble("extMonFNAL.sensor_xoffset", det->m_sensor_xoffset, det->m_sensor_zoffset.size());
      config.getVectorDouble("extMonFNAL.sensor_yoffset", det->m_sensor_yoffset, det->m_sensor_zoffset.size());
      config.getVectorDouble("extMonFNAL.sensor_halfdx",  det->m_sensor_halfdx,  det->m_sensor_zoffset.size());
      config.getVectorDouble("extMonFNAL.sensor_halfdy",  det->m_sensor_halfdy,  det->m_sensor_zoffset.size());
      config.getVectorDouble("extMonFNAL.sensor_halfdz",  det->m_sensor_halfdz,  det->m_sensor_zoffset.size());
      config.getVectorDouble("extMonFNAL.readout_halfdz", det->m_readout_halfdz, det->m_sensor_zoffset.size());

      //----------------------------------------------------------------
      // Test material plates
      const double testMaterialThickness = config.getDouble("extMonFNAL.testMaterial.thickness");
      const double testMaterialXY = config.getDouble("extMonFNAL.testMaterial.halfsize");

      det->m_testMaterialHalfSize.resize(3);
      det->m_testMaterialHalfSize[0] = testMaterialXY;
      det->m_testMaterialHalfSize[1] = testMaterialXY;
      det->m_testMaterialHalfSize[2] = testMaterialThickness/2;

      det->m_testMaterialDistanceToDetector = config.getDouble("extMonFNAL.testMaterial.distanceToDetector");
      det->m_testMaterialPitch = testMaterialThickness + config.getDouble("extMonFNAL.testMaterial.spacing");

      config.getVectorString("extMonFNAL.testMaterial.materials", det->m_testMaterialNames);

      //----------------------------------------------------------------
      // Compute the size of the detector box.  It should contain pixel planes and the test materials

      det->m_detectorHalfSize.resize(3);
      for(unsigned i=0; i<det->m_sensor_zoffset.size(); ++i) {
        det->m_detectorHalfSize[0] =
          std::max(det->m_detectorHalfSize[0], det->m_sensor_halfdx[i] + std::abs(det->m_sensor_xoffset[i]));

        det->m_detectorHalfSize[1] =
          std::max(det->m_detectorHalfSize[1], det->m_sensor_halfdy[i] + std::abs(det->m_sensor_yoffset[i]));

        det->m_detectorHalfSize[2] =
          std::max(det->m_detectorHalfSize[2], det->m_sensor_halfdz[i] + std::abs(det->m_sensor_zoffset[i])
                   + 2*det->m_readout_halfdz[i]);
      }

      //----------------------------------------------------------------
      det->m_detectorRotationInMu2e = room.collimator2RotationInMu2e();

      const double detectorDistanceToWall = config.getDouble("extMonFNAL.detectorDistanceToWall");

      det->m_detectorCenterInMu2e = room.filterExitInMu2e()
        + room.collimator2RotationInMu2e() * CLHEP::Hep3Vector(0, 0, -detectorDistanceToWall);

      //----------------------------------------------------------------
      // Coordinate transform info

      det->m_coordinateCenterInMu2e = room.filterExitInMu2e();
      det->m_coordinateRotationInMu2e = room.collimator2RotationInMu2e().inverse();

      //----------------------------------------------------------------
      if(verbose) {
        std::cout<<"ExtMonFNAL_Maker: detector center in Mu2e = "<<det->m_detectorCenterInMu2e<<std::endl;
        std::cout<<"ExtMonFNAL_Maker: coordinateCenterInMu2e = "<<det->m_coordinateCenterInMu2e<<std::endl;
        std::cout<<"ExtMonFNAL_Maker: coordinateRotationInMu2e = "<<det->m_coordinateRotationInMu2e<<std::endl;
      }

      AGDEBUG("ExtMonFNAL maker end");

      return det;

    } // make()

  } // namespace ExtMonFNAL
} // namespace mu2e
