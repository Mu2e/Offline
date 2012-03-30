// Andrei Gaponenko, 2011

#include "ProtonBeamDumpGeom/inc/ProtonBeamDumpMaker.hh"

#include <algorithm>
#include <iterator>
#include <iostream>
#include <cmath>

#include "cetlib/exception.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include "Mu2eUtilities/inc/SimpleConfig.hh"

#include "ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"

namespace mu2e {
  //================================================================
  ProtonBeamDump::CollimatorExtMonFNAL
  ProtonBeamDumpMaker::readCollimatorExtMonFNAL(const std::string& name,
                                                double angleH,
                                                double angleV,
                                                const SimpleConfig& c) {

    ProtonBeamDump::CollimatorExtMonFNAL col;

    col._name = name;
    col._horizontalLength = c.getDouble("extMonFilter."+name+".horizontalLength");
    c.getVectorDouble("extMonFilter."+name+".channelWidth",            col._channelWidth           , 2);
    c.getVectorDouble("extMonFilter."+name+".channelHeigh",            col._channelHeight          , 2);
    c.getVectorDouble("extMonFilter."+name+".alignmentPlugRadius",     col._alignmentPlugRadius    , 2);
    c.getVectorDouble("extMonFilter."+name+".alignmentHoleRClearance", col._alignmentHoleRClearance, 2);
    col._radiusTransitiondZ = c.getDouble("extMonFilter."+name+".radiusTransitiondZ");
    col._angleH = angleH;
    col._angleV = angleV;

    return col;
  }

  //================================================================
  ProtonBeamDump::FilterMagnetExtMonFNAL
  ProtonBeamDumpMaker::readFilterMagnetExtMonFNAL(const SimpleConfig& c) {
    ProtonBeamDump::FilterMagnetExtMonFNAL mag;
    c.getVectorDouble("extMonFilter.magnet.outerHalfSize", mag._outerHalfSize, 3);
    mag._apertureWidth = c.getDouble("extMonFilter.magnet.apertureWidth") * CLHEP::mm;
    mag._apertureHeight = c.getDouble("extMonFilter.magnet.apertureHeight") * CLHEP::mm;
    mag._fieldStrength = c.getDouble("extMonFilter.magnet.fieldStrength") * CLHEP::tesla;
    return mag;
  }

  //================================================================
  std::auto_ptr<ProtonBeamDump> ProtonBeamDumpMaker::make(const SimpleConfig& c) {

    std::auto_ptr<ProtonBeamDump> m_det(new ProtonBeamDump());

    int verbose = c.getInt("protonBeamDump.verbosityLevel", 0);

    c.getVectorDouble("protonBeamDump.coreHalfSize", m_det->_coreHalfSize, 3);
    c.getVectorDouble("protonBeamDump.neutronCaveHalfSize", m_det->_neutronCaveHalfSize, 3);
    c.getVectorDouble("protonBeamDump.mouthHalfSize", m_det->_mouthHalfSize, 3);
    c.getVectorDouble("extMonFilter.magnetPitHalfSize", m_det->_magnetPitHalfSize, 3);
    m_det->_minCoreShieldingThickness = c.getDouble("protonBeamDump.minCoreShieldingThickness");

    // The default width is determined by core and shielding parameters.   This option lets to increase the default:
    const double enclosureHalfWidthMin = c.getDouble("protonBeamDump.enclosureHalfWidthMin", 0);

    // position
    m_det->_coreCenterInMu2e = c.getHep3Vector("protonBeamDump.coreCenterInMu2e");
    const double coreRotY = m_det->_coreRotY = c.getDouble("protonBeamDump.coreRotY") * CLHEP::degree;

    m_det->_filterEntranceOffsetX = c.getDouble("extMonFilter.entranceOffsetX") * CLHEP::mm;
    m_det->_filterEntranceOffsetY = c.getDouble("extMonFilter.entranceOffsetY") * CLHEP::mm;

    const double angleH = m_det->_filterAngleH =  c.getDouble("extMonFilter.angleH") * CLHEP::radian;
    const double entranceAngleV = m_det->_filterEntranceAngleV =  c.getDouble("extMonFilter.entranceAngleV") * CLHEP::radian;

    const double pNominal = m_det->_extMonFilter_nominalMomentum = c.getDouble("extMonFilter.nominalMomentum") * CLHEP::MeV;

    m_det->_filterMagnet = readFilterMagnetExtMonFNAL(c);


    m_det->_collimator1 = readCollimatorExtMonFNAL("collimator1", angleH, entranceAngleV, c);
    m_det->_collimator2 = readCollimatorExtMonFNAL("collimator2",
                                                   angleH,
                                                   entranceAngleV - 2 * m_det->_filterMagnet.trackBendHalfAngle(pNominal),
                                                   c);

    // Compute the overall size
    m_det->_enclosureHalfSize.resize(3);

    m_det->_enclosureHalfSize[0] =  std::max(m_det->_coreHalfSize[0] + m_det->_minCoreShieldingThickness, enclosureHalfWidthMin);

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
    m_det->_coreCenterInEnclosure[1] = -0.5*m_det->_magnetPitHalfSize[1];
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

    //----------------------------------------------------------------
    // Compute the placement of filter elements

    // collimator1

    const CLHEP::Hep3Vector collimator1CenterInEnclosure(m_det->coreCenterInEnclosure()[0]
                                                         + m_det->filterEntranceOffsetX()
                                                         + 0.5*m_det->collimator1().horizontalLength()*tan(m_det->collimator1().angleH()),

                                                         m_det->coreCenterInEnclosure()[1]
                                                         + m_det->filterEntranceOffsetY()
                                                         + 0.5*m_det->collimator1().horizontalLength()*tan(m_det->collimator1().angleV()),

                                                         m_det->enclosureHalfSize()[2] - 0.5*m_det->collimator1().horizontalLength());

    m_det->_collimator1CenterInEnclosure = collimator1CenterInEnclosure;

    //----------------------------------------------------------------
    // Position the magnet.  Keep B field horizontal.  Then there is
    // no bend in the XZ plane and the magnet center should be on the
    // continuation of the collimator1 center line.
    //
    // The magnet should also have the same rotation in the projection
    // on XZ as the collimator(s).
    //
    // The reference trajectory: a positive particle with the nominal
    // momentum, travelling parallel to collimator1 axis at the bottom
    // center of collimator1 aperture, should enter magnet at the
    // bottom center of the magnet aperture, bend down, and exit a the
    // bottom center of magnet aperture.

    const double magnetAngleV = m_det->_filterMagnetAngleV =
      m_det->filterEntranceAngleV() - m_det->filterMagnet().trackBendHalfAngle(m_det->extMonFilter_nominalMomentum());

    // Arbitrary fix the (center,bottom,center) point of the magnet
    // aperture at the Z position of the magnet room center.

    // Z of magnet entrance at the (center,bottom) of aperture (in enclosure coords)
    const double magnetEntranceZ = m_det->magnetPitCenterInEnclosure()[2]
      + m_det->filterMagnet().outerHalfSize()[2] * cos(magnetAngleV) * cos(m_det->filterAngleH());

    // Height of the entrance point (in enclosure coords)
    const double magnetEntranceY = m_det->coreCenterInEnclosure()[1] + m_det->filterEntranceOffsetY()
      + tan(m_det->filterEntranceAngleV()) * (m_det->enclosureHalfSize()[2] - magnetEntranceZ)

      // the entrance point should be on the trajectory of reference
      // particle at the *bottom* of the collimator1 channel - shift for that:
      - 0.5 * m_det->collimator1().channelHeight()[m_det->collimator1().channelHeight().size()-1]/cos(m_det->collimator1().angleV())
      ;

    // X of the entrance point (in enclosure coords)
    const double magnetEntranceX = m_det->coreCenterInEnclosure()[0] + m_det->filterEntranceOffsetX()
      + tan(m_det->filterAngleH()) * (m_det->enclosureHalfSize()[2] - magnetEntranceZ);

    // length from magnet entrance to magnet center in projection on (XZ)
    const double projLength =
      (m_det->filterMagnet().outerHalfSize()[2]* cos(magnetAngleV) - 0.5*m_det->filterMagnet().apertureHeight()*sin(magnetAngleV));

    // Compute the center of the magnet
    const CLHEP::Hep3Vector magnetCenterInEnclosure
      (magnetEntranceX + sin(m_det->filterAngleH()) * projLength
       ,
       magnetEntranceY + sin(magnetAngleV)*m_det->filterMagnet().outerHalfSize()[2] + cos(magnetAngleV)*0.5*m_det->filterMagnet().apertureHeight()
       ,
       magnetEntranceZ - cos(m_det->filterAngleH()) * projLength
       );

    m_det->_filterMagnetCenterInEnclosure = magnetCenterInEnclosure;

    //----------------------------------------------------------------
    // collimator2: position so that the reference trajectory enters it at the bottom center
    // of the aperture

    const double col2CenterZ = -m_det->_enclosureHalfSize[2] + 0.5*m_det->_collimator2.horizontalLength();

    const double col2CenterX = m_det->_coreCenterInEnclosure[0] + m_det->_filterEntranceOffsetX +
      tan(m_det->filterAngleH()) * (m_det->_enclosureHalfSize[2] - col2CenterZ);

    // Z of the exit point from the magnet, at the (center,bottom) of aperture (in enclosure coords)
    const double magnetExitZ = m_det->magnetPitCenterInEnclosure()[2]
      - m_det->filterMagnet().outerHalfSize()[2] * cos(magnetAngleV) * cos(m_det->filterAngleH());

    // Height of the exit point (in enclosure coords)
    const double magnetExitY = magnetEntranceY + sin(magnetAngleV) * 2*m_det->_filterMagnet._outerHalfSize[2];

    const double col2CenterY = magnetExitY + tan(m_det->_collimator2.angleV())*(magnetExitZ - col2CenterZ)
      + 0.5 * m_det->_collimator2._channelHeight[0]/cos(m_det->_collimator2.angleV());

    m_det->_collimator2CenterInEnclosure = CLHEP::Hep3Vector(col2CenterX, col2CenterY, col2CenterZ);

    //----------------------------------------------------------------
    if(verbose) {
      std::cout<<__func__<<": coreOffset = "<<coreOffset<<std::endl;
      std::cout<<__func__<<": ProtonBeamDump core center in mu2e = "<<m_det->_coreCenterInMu2e<<std::endl;
      std::cout<<__func__<<": ProtonBeamDump enclosure center in mu2e = "<<m_det->_enclosureCenterInMu2e<<std::endl;
      std::cout<<__func__<<": coreCenterInEnclosure = "<<m_det->_coreCenterInEnclosure<<std::endl;
      std::cout<<__func__<<": enclosure half size = { "
               <<m_det->_enclosureHalfSize[0]<<", "
               <<m_det->_enclosureHalfSize[1]<<", "
               <<m_det->_enclosureHalfSize[2]<<" }"
               <<std::endl;
      std::cout<<__func__<<": magnetPitCenterInEnclosure = "<<m_det->_magnetPitCenterInEnclosure<<std::endl;
      std::cout<<__func__<<": shieldingFaceXmin = "<<m_det->_shieldingFaceXmin
               <<", Xmax = "<<m_det->_shieldingFaceXmax<<std::endl;
      std::cout<<__func__<<": shieldingFaceZatXmin = "<<m_det->_shieldingFaceZatXmin
               <<", ZatXmax = "<<m_det->_shieldingFaceZatXmax<<std::endl;
      std::cout<<__func__<<": filter nominal momentum = "<<m_det->extMonFilter_nominalMomentum()/CLHEP::GeV<<" GeV/c"<<std::endl;
      std::cout<<__func__<<": filterAngleH = "<<m_det->filterAngleH()<<std::endl;
      std::cout<<__func__<<": filter half bend angle  = "<<m_det->filterMagnet().trackBendHalfAngle(m_det->extMonFilter_nominalMomentum())<<std::endl;
      std::cout<<__func__<<": filter.angleV = "<<m_det->filterEntranceAngleV()
               <<", c1.angleV  = "<<m_det->collimator1().angleV()
               <<", magnet.angleV = "<<m_det->filterMagnetAngleV()
               <<", c2.angleV() = "<<m_det->collimator2().angleV()<<std::endl;
      std::cout<<__func__<<": collimator1CenterInEnclosure = "<<m_det->_collimator1CenterInEnclosure<<std::endl;
      std::cout<<__func__<<": collimator1.horizontalLength = "<<m_det->_collimator1._horizontalLength<<std::endl;
      std::cout<<__func__<<": collimator2CenterInEnclosure = "<<m_det->_collimator2CenterInEnclosure<<std::endl;
      std::cout<<__func__<<": collimator2.horizontalLength = "<<m_det->_collimator2._horizontalLength<<std::endl;
    }

    return m_det;
  } // make()

} // namespace mu2e
