// Andrei Gaponenko, 2011

#include "GeometryService/inc/ProtonBeamDumpMaker.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/ProductionTarget.hh"

#include <algorithm>
#include <iterator>
#include <iostream>
#include <cmath>

#include "cetlib/exception.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include "Mu2eUtilities/inc/SimpleConfig.hh"

#include "GeometryService/inc/ProtonBeamDump.hh"

namespace mu2e {
  
  ProtonBeamDumpMaker::ProtonBeamDumpMaker(const SimpleConfig& c) 
    : m_det(new ProtonBeamDump())
  {
    int verbose = c.getInt("protonBeamDump.verbosityLevel", 0);
    
    c.getVectorDouble("protonBeamDump.coreHalfSize", m_det->_coreHalfSize, 3);
    c.getVectorDouble("protonBeamDump.neutronCaveHalfSize", m_det->_neutronCaveHalfSize, 3);
    c.getVectorDouble("protonBeamDump.mouthHalfSize", m_det->_mouthHalfSize, 3);
    c.getVectorDouble("extMonFilter.magnetPitHalfSize", m_det->_magnetPitHalfSize, 3);
    m_det->_minCoreShieldingThickness = c.getDouble("protonBeamDump.minCoreShieldingThickness");

    m_det->_collimator1horizLength = c.getDouble("extMonFilter.collimator1.horizontalLength");
    m_det->_collimator2horizLength = c.getDouble("extMonFilter.collimator2.horizontalLength");

    // position
    m_det->_coreCenterInMu2e = c.getHep3Vector("protonBeamDump.coreCenterInMu2e");
    m_det->_coreRotY = c.getDouble("protonBeamDump.coreRotY") * CLHEP::degree;

    // Compute the overall size
    m_det->_enclosureHalfSize.resize(3);

    m_det->_enclosureHalfSize[0] =  m_det->_coreHalfSize[0] + m_det->_minCoreShieldingThickness;

    m_det->_enclosureHalfSize[1] =  m_det->_coreHalfSize[1] + m_det->_minCoreShieldingThickness + m_det->_magnetPitHalfSize[1];

    const double dumpPlusShieldingHalfLength =
      m_det->_mouthHalfSize[2] + m_det->_neutronCaveHalfSize[2] + m_det->_coreHalfSize[2] + 0.5*m_det->_minCoreShieldingThickness;
    
    const double totalFilterLength = 
      m_det->_collimator1horizLength + 2*m_det->_magnetPitHalfSize[2] + m_det->_collimator2horizLength;
    
    m_det->_enclosureHalfSize[2] = std::max(dumpPlusShieldingHalfLength, 0.5*totalFilterLength);

    if(verbose) {
      std::cout<<__func__<<": ProtonBeamDump enclosure half size = ";
      std::copy(m_det->_enclosureHalfSize.begin(), m_det->_enclosureHalfSize.end(), std::ostream_iterator<double>(std::cout, ", "));
      std::cout<<std::endl;
    }

    // compute the position of the overall enclosure
    m_det->_enclosureRotationInMu2e.rotateY(-m_det->_coreRotY);

    // The offset of the enclosure center w.r.t. the core, along the dump z
    const double coreOffset = 2*m_det->_mouthHalfSize[2] + 2*m_det->_neutronCaveHalfSize[2]
      + m_det->_coreHalfSize[2] - m_det->_enclosureHalfSize[2];

    m_det->_enclosureCenterInMu2e[0] = m_det->_coreCenterInMu2e[0] + coreOffset*sin(m_det->_coreRotY);
    m_det->_enclosureCenterInMu2e[1] = m_det->_coreCenterInMu2e[1] + 0.5*m_det->_magnetPitHalfSize[1];
    m_det->_enclosureCenterInMu2e[2] = m_det->_coreCenterInMu2e[2] + coreOffset*cos(m_det->_coreRotY);

    // core relative to the enclosure
    m_det->_coreCenterInEnclosure[0] = 0;
    m_det->_coreCenterInEnclosure[1] = -m_det->_magnetPitHalfSize[1];
    m_det->_coreCenterInEnclosure[2] = -coreOffset;

    // position of the magnet pit
    m_det->_magnetPitCenterInEnclosure[0] = 0.;
    m_det->_magnetPitCenterInEnclosure[1] = m_det->_enclosureHalfSize[1] - m_det->_magnetPitHalfSize[1];
    m_det->_magnetPitCenterInEnclosure[2] = m_det->_enclosureHalfSize[2] - m_det->_collimator1horizLength - m_det->_magnetPitHalfSize[2];

    if(verbose) {
      std::cout<<__func__<<": coreOffset = "<<coreOffset<<std::endl;
      std::cout<<__func__<<": ProtonBeamDump core center in mu2e = "<<m_det->_coreCenterInMu2e<<std::endl;
      std::cout<<__func__<<": ProtonBeamDump enclosure center in mu2e = "<<m_det->_enclosureCenterInMu2e<<std::endl;
      std::cout<<__func__<<": coreCenterInEnclosure = "<<m_det->_coreCenterInEnclosure<<std::endl;
      std::cout<<__func__<<": magnetPitCenterInEnclosure = "<<m_det->_magnetPitCenterInEnclosure<<std::endl;
    }
  }
}
