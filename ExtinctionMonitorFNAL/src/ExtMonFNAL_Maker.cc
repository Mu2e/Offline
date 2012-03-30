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

      std::auto_ptr<ExtMon> m_det(new ExtMon(hs, config.getDouble("extmon_fnal.roomCenterHeightAboveDumpCore")));

      m_det->m_roomCenterInMu2e[0] = dump.enclosureCenterInMu2e()[0]
        - sin(dump.coreRotY())*(dump.enclosureHalfSize()[2]+m_det->m_roomHalfSize[2]);

      m_det->m_roomCenterInMu2e[1] = dump.coreCenterInMu2e()[1] + m_det->m_roomCenterHeightAboveDumpCore;

      m_det->m_roomCenterInMu2e[2] = dump.enclosureCenterInMu2e()[2]
        - cos(dump.coreRotY())*(dump.enclosureHalfSize()[2]+m_det->m_roomHalfSize[2]);


      config.getVectorDouble("extmon_fnal.sensor_zoffset", m_det->m_sensor_zoffset, -1);
      config.getVectorDouble("extmon_fnal.sensor_xoffset", m_det->m_sensor_xoffset, m_det->m_sensor_zoffset.size());
      config.getVectorDouble("extmon_fnal.sensor_yoffset", m_det->m_sensor_yoffset, m_det->m_sensor_zoffset.size());
      config.getVectorDouble("extmon_fnal.sensor_halfdx",  m_det->m_sensor_halfdx,  m_det->m_sensor_zoffset.size());
      config.getVectorDouble("extmon_fnal.sensor_halfdy",  m_det->m_sensor_halfdy,  m_det->m_sensor_zoffset.size());
      config.getVectorDouble("extmon_fnal.sensor_halfdz",  m_det->m_sensor_halfdz,  m_det->m_sensor_zoffset.size());
      config.getVectorDouble("extmon_fnal.readout_halfdz", m_det->m_readout_halfdz, m_det->m_sensor_zoffset.size());

      m_det->m_detectorHalfSize.resize(3);
      for(unsigned i=0; i<m_det->m_sensor_zoffset.size(); ++i) {
        m_det->m_detectorHalfSize[0] =
          std::max(m_det->m_detectorHalfSize[0], m_det->m_sensor_halfdx[i] + std::abs(m_det->m_sensor_xoffset[i]));

        m_det->m_detectorHalfSize[1] =
          std::max(m_det->m_detectorHalfSize[1], m_det->m_sensor_halfdy[i] + std::abs(m_det->m_sensor_yoffset[i]));

        m_det->m_detectorHalfSize[2] =
          std::max(m_det->m_detectorHalfSize[2], m_det->m_sensor_halfdz[i] + std::abs(m_det->m_sensor_zoffset[i])
                   + 2*m_det->m_readout_halfdz[i]);


      }

      //----------------------------------------------------------------
      m_det->m_detectorRotationInRoom
        .rotateX(-dump.collimator2().angleV())
        .rotateY(+dump.collimator2().angleH());

      const double detectorDistanceToWall = config.getDouble("extmon_fnal.detectorDistanceToWall");
      const double detectorCenterInRoomZ = m_det->m_roomHalfSize[2]
        - (detectorDistanceToWall + m_det->m_detectorHalfSize[2]) * cos(dump.collimator2().angleV()) * cos(dump.collimator2().angleH());

      const double col2CenterInRoomZ = dump.collimator2CenterInEnclosure()[2] + m_det->m_roomHalfSize[2] + dump.enclosureHalfSize()[2];
      const double col2CenterInRoomX = dump.collimator2CenterInEnclosure()[0];
      const double col2CenterInRoomY = dump.collimator2CenterInEnclosure()[1] + (dump.enclosureCenterInMu2e()[1] - m_det->m_roomCenterInMu2e[1]);

      const double detectorCenterInRoomX = col2CenterInRoomX - tan(dump.collimator2().angleH()) * (detectorCenterInRoomZ - col2CenterInRoomZ);
      const double detectorCenterInRoomY = col2CenterInRoomY - tan(dump.collimator2().angleV()) * (detectorCenterInRoomZ - col2CenterInRoomZ);

      m_det->m_detectorCenterInRoom = CLHEP::Hep3Vector(detectorCenterInRoomX, detectorCenterInRoomY, detectorCenterInRoomZ);

      //----------------------------------------------------------------
      if(verbose) {
        std::cout<<__func__<<": ExtMonFNAL room center in mu2e = "<<m_det->roomCenterInMu2e()<<std::endl;
        std::cout<<__func__<<": ExtMonFNAL room rotY = "<<dump.coreRotY()<<std::endl;
        std::cout<<__func__<<": ExtMonFNAL room half size = ("
                 <<m_det->roomHalfSize()[0]<<", "
                 <<m_det->roomHalfSize()[1]<<", "
                 <<m_det->roomHalfSize()[2]<<")"
                 <<std::endl;

        std::cout<<__func__<<": ExtMonFNAL detector center in room = "<<m_det->m_detectorCenterInRoom<<std::endl;
      }

      return m_det;

    } // make()

  } // namespace ExtMonFNAL
} // namespace mu2e
