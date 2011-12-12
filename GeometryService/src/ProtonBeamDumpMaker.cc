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
  //================================================================
  ProtonBeamDump::CollimatorExtMonFNAL 
  ProtonBeamDumpMaker::readCollimatorExtMonFNAL(const std::string& name, const SimpleConfig& c) {

    ProtonBeamDump::CollimatorExtMonFNAL col;
    
    col._name = name;
    col._horizontalLength = c.getDouble("extMonFilter."+name+".horizontalLength");
    c.getVectorDouble("extMonFilter."+name+".channelWidth",            col._channelWidth           , 2);
    c.getVectorDouble("extMonFilter."+name+".channelHeigh",            col._channelHeight          , 2);
    c.getVectorDouble("extMonFilter."+name+".alignmentPlugRadius",     col._alignmentPlugRadius    , 2);
    c.getVectorDouble("extMonFilter."+name+".alignmentHoleRClearance", col._alignmentHoleRClearance, 2);
    col._radiusTransitiondZ = c.getDouble("extMonFilter."+name+".radiusTransitiondZ");
    col._angleH = c.getDouble("extMonFilter."+name+".angleH");
    col._angleV = c.getDouble("extMonFilter."+name+".angleV");
    col._entranceOffsetX = c.getDouble("extMonFilter."+name+".entranceOffsetX");
    col._entranceOffsetY = c.getDouble("extMonFilter."+name+".entranceOffsetY");

    return col;
  }

  //================================================================
  ProtonBeamDumpMaker::ProtonBeamDumpMaker(const SimpleConfig& c) 
    : m_det(new ProtonBeamDump())
  {
    int verbose = c.getInt("protonBeamDump.verbosityLevel", 0);
    
    c.getVectorDouble("protonBeamDump.coreHalfSize", m_det->_coreHalfSize, 3);
    c.getVectorDouble("protonBeamDump.neutronCaveHalfSize", m_det->_neutronCaveHalfSize, 3);
    c.getVectorDouble("protonBeamDump.mouthHalfSize", m_det->_mouthHalfSize, 3);
    c.getVectorDouble("extMonFilter.magnetPitHalfSize", m_det->_magnetPitHalfSize, 3);
    m_det->_minCoreShieldingThickness = c.getDouble("protonBeamDump.minCoreShieldingThickness");

    m_det->_collimator1 = readCollimatorExtMonFNAL("collimator1", c);
    m_det->_collimator2 = readCollimatorExtMonFNAL("collimator2", c);

    // position
    m_det->_coreCenterInMu2e = c.getHep3Vector("protonBeamDump.coreCenterInMu2e");
    const double coreRotY = m_det->_coreRotY = c.getDouble("protonBeamDump.coreRotY") * CLHEP::degree;

    // Compute the overall size
    m_det->_enclosureHalfSize.resize(3);

    m_det->_enclosureHalfSize[0] =  m_det->_coreHalfSize[0] + m_det->_minCoreShieldingThickness;

    m_det->_enclosureHalfSize[1] =  m_det->_coreHalfSize[1] + m_det->_minCoreShieldingThickness + m_det->_magnetPitHalfSize[1];

    const double dumpPlusShieldingHalfLength =
      m_det->_mouthHalfSize[2] + m_det->_neutronCaveHalfSize[2] + m_det->_coreHalfSize[2] + 0.5*m_det->_minCoreShieldingThickness;
    
    const double totalFilterLength = 
      m_det->_collimator1._horizontalLength + 2*m_det->_magnetPitHalfSize[2] + m_det->_collimator2._horizontalLength;
    
    m_det->_enclosureHalfSize[2] = std::max(dumpPlusShieldingHalfLength, 0.5*totalFilterLength);

    // shorthand notation
    const std::vector<double>& enclosureHalfSize = m_det->_enclosureHalfSize;

    if(verbose) {
      std::cout<<__func__<<": ProtonBeamDump enclosure half size = ";
      std::copy(enclosureHalfSize.begin(), enclosureHalfSize.end(), std::ostream_iterator<double>(std::cout, ", "));
      std::cout<<std::endl;
    }

    // compute the position of the overall enclosure
    m_det->_enclosureRotationInMu2e.rotateY(-coreRotY);

    // The offset of the enclosure center w.r.t. the core, along the dump z
    const double coreOffset = 2*m_det->_mouthHalfSize[2] + 2*m_det->_neutronCaveHalfSize[2]
      + m_det->_coreHalfSize[2] - enclosureHalfSize[2];

    m_det->_enclosureCenterInMu2e[0] = m_det->_coreCenterInMu2e[0] + coreOffset*sin(coreRotY);
    m_det->_enclosureCenterInMu2e[1] = m_det->_coreCenterInMu2e[1] + 0.5*m_det->_magnetPitHalfSize[1];
    m_det->_enclosureCenterInMu2e[2] = m_det->_coreCenterInMu2e[2] + coreOffset*cos(coreRotY);

    const CLHEP::Hep3Vector& enclosureCenterInMu2e = m_det->_enclosureCenterInMu2e;

    // core relative to the enclosure
    m_det->_coreCenterInEnclosure[0] = 0;
    m_det->_coreCenterInEnclosure[1] = -m_det->_magnetPitHalfSize[1];
    m_det->_coreCenterInEnclosure[2] = -coreOffset;

    // position of the magnet pit
    m_det->_magnetPitCenterInEnclosure[0] = 0.;
    m_det->_magnetPitCenterInEnclosure[1] = enclosureHalfSize[1] - m_det->_magnetPitHalfSize[1];
    m_det->_magnetPitCenterInEnclosure[2] = enclosureHalfSize[2] - m_det->_collimator1._horizontalLength - m_det->_magnetPitHalfSize[2];

    // Shielding face coordinates
    m_det->_shieldingFaceXmin = enclosureCenterInMu2e[0] 
      + enclosureHalfSize[2] * sin(coreRotY)
      - enclosureHalfSize[0] * cos(coreRotY)
      ;

    m_det->_shieldingFaceXmax = enclosureCenterInMu2e[0] 
      + enclosureHalfSize[2] * sin(coreRotY)
      + enclosureHalfSize[0] * cos(coreRotY)
      ;

    m_det->_shieldingFaceZatXmin = enclosureCenterInMu2e[2]
      + enclosureHalfSize[2] * cos(coreRotY)
      + enclosureHalfSize[0] * sin(coreRotY)
      ;

    m_det->_shieldingFaceZatXmax = enclosureCenterInMu2e[2]
      + enclosureHalfSize[2] * cos(coreRotY)
      - enclosureHalfSize[0] * sin(coreRotY)
      ;

    if(verbose) {
      std::cout<<__func__<<": coreOffset = "<<coreOffset<<std::endl;
      std::cout<<__func__<<": ProtonBeamDump core center in mu2e = "<<m_det->_coreCenterInMu2e<<std::endl;
      std::cout<<__func__<<": ProtonBeamDump enclosure center in mu2e = "<<m_det->_enclosureCenterInMu2e<<std::endl;
      std::cout<<__func__<<": coreCenterInEnclosure = "<<m_det->_coreCenterInEnclosure<<std::endl;
      std::cout<<__func__<<": magnetPitCenterInEnclosure = "<<m_det->_magnetPitCenterInEnclosure<<std::endl;
    }
  }
}
