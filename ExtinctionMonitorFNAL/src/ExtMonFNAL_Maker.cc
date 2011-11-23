// Andrei Gaponenko, 2011

#include "ExtinctionMonitorFNAL/inc/ExtMonFNAL_Maker.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/ProtonBeamDump.hh"

#include "cetlib/exception.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "ExtinctionMonitorFNAL/inc/ExtMonFNAL.hh"

namespace mu2e {
  namespace ExtMonFNAL {

    ExtMonMaker::ExtMonMaker(const SimpleConfig& config) 
      : m_det(0)
    {
      std::vector<double> hs;
      config.getVectorDouble("extmon_fnal.roomHalfSize", hs, 3);	
      m_det.reset(new ExtMon(hs, config.getDouble("extmon_fnal.roomCenterHeightAboveDumpCore")));

      GeomHandle<ProtonBeamDump> dump;
      m_det->m_roomCenterInMu2e[0] = dump->enclosureCenterInMu2e()[0]
	- sin(dump->coreRotY())*(dump->enclosureHalfSize()[2]+m_det->m_roomHalfSize[2]);

      m_det->m_roomCenterInMu2e[1] = dump->coreCenterInMu2e()[1] + m_det->m_roomCenterHeightAboveDumpCore;

      m_det->m_roomCenterInMu2e[2] = dump->enclosureCenterInMu2e()[2]
	- cos(dump->coreRotY())*(dump->enclosureHalfSize()[2]+m_det->m_roomHalfSize[2]);

      m_det->m_detectorCenterInRoom = config.getHep3Vector("extmon_fnal.detectorCenterInRoom");
      
      config.getVectorDouble("extmon_fnal.sensor_zoffset", m_det->m_sensor_zoffset, -1);
      config.getVectorDouble("extmon_fnal.sensor_xoffset", m_det->m_sensor_xoffset, m_det->m_sensor_zoffset.size());
      config.getVectorDouble("extmon_fnal.sensor_yoffset", m_det->m_sensor_yoffset, m_det->m_sensor_zoffset.size());
      config.getVectorDouble("extmon_fnal.sensor_halfdx",  m_det->m_sensor_halfdx,  m_det->m_sensor_zoffset.size());
      config.getVectorDouble("extmon_fnal.sensor_halfdy",  m_det->m_sensor_halfdy,  m_det->m_sensor_zoffset.size());
      config.getVectorDouble("extmon_fnal.sensor_halfdz",  m_det->m_sensor_halfdz,  m_det->m_sensor_zoffset.size());
      config.getVectorDouble("extmon_fnal.readout_halfdz", m_det->m_readout_halfdz, m_det->m_sensor_zoffset.size());

      m_det->m_detectorRotationInRoom
	.rotateX(config.getDouble("extmon_fnal.detectorRotationInRoomX")*CLHEP::degree)
	.rotateY(config.getDouble("extmon_fnal.detectorRotationInRoomY")*CLHEP::degree);
      
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
      
    }

  }
}
