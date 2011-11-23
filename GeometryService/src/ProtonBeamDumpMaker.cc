// Andrei Gaponenko, 2011

#include "GeometryService/inc/ProtonBeamDumpMaker.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/ProductionTarget.hh"
#include "GeometryService/inc/Mu2eBuilding.hh"

#include <algorithm>
#include <cmath>

#include "cetlib/exception.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include "Mu2eUtilities/inc/SimpleConfig.hh"

#include "GeometryService/inc/ProtonBeamDump.hh"

namespace mu2e {
  
  ProtonBeamDumpMaker::ProtonBeamDumpMaker(const SimpleConfig& c) 
    : m_det(new ProtonBeamDump())
  {
    
    c.getVectorDouble("protonBeamDump.coreHalfSize", m_det->_coreHalfSize, 3);
    c.getVectorDouble("protonBeamDump.neutronCaveHalfSize", m_det->_neutronCaveHalfSize, 3);
    c.getVectorDouble("protonBeamDump.mouthHalfSize", m_det->_mouthHalfSize, 3);
    c.getVectorDouble("extMonFilter.magnetHollowHalfSize", m_det->_magnetHollowHalfSize, 3);
    m_det->_minCoreShieldingThickness = c.getDouble("protonBeamDump.minCoreShieldingThickness");

    m_det->_collimator1horizLength = c.getDouble("extMonFilter.collimator1.horizontalLength");
    m_det->_collimator2horizLength = c.getDouble("extMonFilter.collimator2.horizontalLength");

    // position
    m_det->_coreCenterInMu2e = c.getHep3Vector("protonBeamDump.coreCenterInMu2e");
    m_det->_coreRotY = c.getDouble("protonBeamDump.coreRotY") * CLHEP::degree;

    // Compute the overall size
    m_det->_enclosureHalfSize.resize(3);

    m_det->_enclosureHalfSize[0] =  m_det->_coreHalfSize[0] + m_det->_minCoreShieldingThickness;

    m_det->_enclosureHalfSize[1] =  m_det->_coreHalfSize[1] + m_det->_minCoreShieldingThickness + m_det->_magnetHollowHalfSize[1];

    const double dumpPlusShieldingLength =
      m_det->_mouthHalfSize[2] + m_det->_neutronCaveHalfSize[2] + m_det->_coreHalfSize[2] + 0.5*m_det->_minCoreShieldingThickness;
    
    const double totalFilterLength = 
      0.5*m_det->_collimator1horizLength + m_det->_magnetHollowHalfSize[2] + 0.5*m_det->_collimator2horizLength;
    
    m_det->_enclosureHalfSize[2] = std::max(dumpPlusShieldingLength, totalFilterLength);

    // compute the position
    m_det->_enclosureRotationInMu2e.rotateY(-m_det->_coreRotY);

    m_det->_enclosureCenterInMu2e[0] = m_det->_coreCenterInMu2e[0];
    m_det->_enclosureCenterInMu2e[1] = m_det->_coreCenterInMu2e[1] + 0.5*m_det->_magnetHollowHalfSize[1];
    m_det->_enclosureCenterInMu2e[2] = m_det->_coreCenterInMu2e[2] - std::max(0., 0.5*(totalFilterLength - dumpPlusShieldingLength));
  }
}
