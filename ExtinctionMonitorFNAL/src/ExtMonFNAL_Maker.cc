#include "ExtinctionMonitorFNAL/inc/ExtMonFNAL_Maker.hh"

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
      config.getVectorDouble("extmon_fnal.logicalEnclosureHalfDim", hs, 3);	
      m_det.reset(new ExtMon(hs, config.getHep3Vector("extmon_fnal.offsetInParent")));
      
      config.getVectorDouble("extmon_fnal.sensor_zoffset", m_det->m_sensor_zoffset, -1);
      config.getVectorDouble("extmon_fnal.sensor_xoffset", m_det->m_sensor_xoffset, m_det->m_sensor_zoffset.size());
      config.getVectorDouble("extmon_fnal.sensor_yoffset", m_det->m_sensor_yoffset, m_det->m_sensor_zoffset.size());
      config.getVectorDouble("extmon_fnal.sensor_halfdx",  m_det->m_sensor_halfdx,  m_det->m_sensor_zoffset.size());
      config.getVectorDouble("extmon_fnal.sensor_halfdy",  m_det->m_sensor_halfdy,  m_det->m_sensor_zoffset.size());
      config.getVectorDouble("extmon_fnal.sensor_halfdz",  m_det->m_sensor_halfdz,  m_det->m_sensor_zoffset.size());
      config.getVectorDouble("extmon_fnal.readout_halfdz", m_det->m_readout_halfdz, m_det->m_sensor_zoffset.size());

      m_det->m_rotationInParent
	.rotateX(config.getDouble("extmon_fnal.rotationX")*CLHEP::degree)
	.rotateY(config.getDouble("extmon_fnal.rotationY")*CLHEP::degree);
    }

  }
}
