// Andrei Gaponenko, 2011

#include "ExtinctionMonitorFNAL/inc/ExtMonFNAL_Maker.hh"

#include "ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"

#include "cetlib/exception.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "ExtinctionMonitorFNAL/inc/ExtMonFNAL.hh"

namespace mu2e {
  namespace ExtMonFNAL {

    std::auto_ptr<ExtMon> ExtMonMaker::make(const SimpleConfig& config, const ProtonBeamDump& dump)
    {
      const int verbose = config.getInt("extmon_fnal.verbosityLevel", 0);

      std::vector<double> hs;
      config.getVectorDouble("extmon_fnal.roomHalfSize", hs, 3);

      std::auto_ptr<ExtMon> det(new ExtMon(hs, config.getDouble("extmon_fnal.roomCenterHeightAboveDumpCore")));

      det->m_roomCenterInMu2e[0] = dump.enclosureCenterInMu2e()[0]
        - sin(dump.coreRotY())*(dump.enclosureHalfSize()[2]+det->m_roomHalfSize[2]);

      det->m_roomCenterInMu2e[1] = dump.coreCenterInMu2e()[1] + det->m_roomCenterHeightAboveDumpCore;

      det->m_roomCenterInMu2e[2] = dump.enclosureCenterInMu2e()[2]
        - cos(dump.coreRotY())*(dump.enclosureHalfSize()[2]+det->m_roomHalfSize[2]);


      config.getVectorDouble("extmon_fnal.sensor_zoffset", det->m_sensor_zoffset, -1);
      config.getVectorDouble("extmon_fnal.sensor_xoffset", det->m_sensor_xoffset, det->m_sensor_zoffset.size());
      config.getVectorDouble("extmon_fnal.sensor_yoffset", det->m_sensor_yoffset, det->m_sensor_zoffset.size());
      config.getVectorDouble("extmon_fnal.sensor_halfdx",  det->m_sensor_halfdx,  det->m_sensor_zoffset.size());
      config.getVectorDouble("extmon_fnal.sensor_halfdy",  det->m_sensor_halfdy,  det->m_sensor_zoffset.size());
      config.getVectorDouble("extmon_fnal.sensor_halfdz",  det->m_sensor_halfdz,  det->m_sensor_zoffset.size());
      config.getVectorDouble("extmon_fnal.readout_halfdz", det->m_readout_halfdz, det->m_sensor_zoffset.size());

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
      det->m_detectorRotationInRoom
        .rotateX(-dump.collimator2().angleV())
        .rotateY(+dump.collimator2().angleH());

      const double detectorDistanceToWall = config.getDouble("extmon_fnal.detectorDistanceToWall");
      const double detectorCenterInRoomZ = det->m_roomHalfSize[2]
        - (detectorDistanceToWall + det->m_detectorHalfSize[2]) * cos(dump.collimator2().angleV()) * cos(dump.collimator2().angleH());

      const double col2CenterInRoomZ = dump.collimator2CenterInEnclosure()[2] + det->m_roomHalfSize[2] + dump.enclosureHalfSize()[2];
      const double col2CenterInRoomX = dump.collimator2CenterInEnclosure()[0];
      const double col2CenterInRoomY = dump.collimator2CenterInEnclosure()[1] + (dump.enclosureCenterInMu2e()[1] - det->m_roomCenterInMu2e[1]);

      const double detectorCenterInRoomX = col2CenterInRoomX - tan(dump.collimator2().angleH()) * (detectorCenterInRoomZ - col2CenterInRoomZ);
      const double detectorCenterInRoomY = col2CenterInRoomY - tan(dump.collimator2().angleV()) * (detectorCenterInRoomZ - col2CenterInRoomZ);

      det->m_detectorCenterInRoom = CLHEP::Hep3Vector(detectorCenterInRoomX, detectorCenterInRoomY, detectorCenterInRoomZ);

      //----------------------------------------------------------------
      if(verbose) {
        std::cout<<__func__<<": ExtMonFNAL room center in mu2e = "<<det->roomCenterInMu2e()<<std::endl;
        std::cout<<__func__<<": ExtMonFNAL room rotY = "<<dump.coreRotY()<<std::endl;
        std::cout<<__func__<<": ExtMonFNAL room half size = ("
                 <<det->roomHalfSize()[0]<<", "
                 <<det->roomHalfSize()[1]<<", "
                 <<det->roomHalfSize()[2]<<")"
                 <<std::endl;

        std::cout<<__func__<<": ExtMonFNAL detector center in room = "<<det->m_detectorCenterInRoom<<std::endl;
      }

      return det;

    } // make()

  } // namespace ExtMonFNAL
} // namespace mu2e
