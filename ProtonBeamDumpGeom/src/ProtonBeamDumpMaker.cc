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
  namespace {
    CLHEP::Hep3Vector enclosureToMu2ePoint(const ProtonBeamDump& dump, double xe, double ye, double ze) {
      const CLHEP::Hep3Vector& ec(dump.enclosureCenterInMu2e());
      return CLHEP::Hep3Vector(  ec[0] + xe*cos(dump.coreRotY()) + ze*sin(dump.coreRotY())
                               , ec[1] + ye
                               , ec[2] - xe*sin(dump.coreRotY()) + ze*cos(dump.coreRotY())
                               );
    }
  }

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

    std::auto_ptr<ProtonBeamDump> dump(new ProtonBeamDump());

    int verbose = c.getInt("protonBeamDump.verbosityLevel", 0);

    c.getVectorDouble("protonBeamDump.coreHalfSize", dump->_coreHalfSize, 3);
    c.getVectorDouble("protonBeamDump.neutronCaveHalfSize", dump->_neutronCaveHalfSize, 3);
    c.getVectorDouble("protonBeamDump.mouthHalfSize", dump->_mouthHalfSize, 3);
    c.getVectorDouble("extMonFilter.magnetPitHalfSize", dump->_magnetPitHalfSize, 3);
    dump->_minCoreShieldingThickness = c.getDouble("protonBeamDump.minCoreShieldingThickness");

    // The default width is determined by core and shielding parameters.   This option lets to increase the default:
    const double enclosureHalfWidthMin = c.getDouble("protonBeamDump.enclosureHalfWidthMin", 0);

    // position
    dump->_coreCenterInMu2e = c.getHep3Vector("protonBeamDump.coreCenterInMu2e");
    const double coreRotY = dump->_coreRotY = c.getDouble("protonBeamDump.coreRotY") * CLHEP::degree;

    dump->_filterEntranceOffsetX = c.getDouble("extMonFilter.entranceOffsetX") * CLHEP::mm;
    dump->_filterEntranceOffsetY = c.getDouble("extMonFilter.entranceOffsetY") * CLHEP::mm;

    const double angleH = dump->_filterAngleH =  c.getDouble("extMonFilter.angleH") * CLHEP::radian;
    const double entranceAngleV = dump->_filterEntranceAngleV =  c.getDouble("extMonFilter.entranceAngleV") * CLHEP::radian;

    const double pNominal = dump->_extMonFilter_nominalMomentum = c.getDouble("extMonFilter.nominalMomentum") * CLHEP::MeV;

    dump->_filterMagnet = readFilterMagnetExtMonFNAL(c);

    dump->_collimator1 = readCollimatorExtMonFNAL("collimator1", angleH, entranceAngleV, c);
    dump->_collimator2 = readCollimatorExtMonFNAL("collimator2",
                                                   angleH,
                                                   entranceAngleV - 2 * dump->_filterMagnet.trackBendHalfAngle(pNominal),
                                                   c);

    // Compute the overall size
    dump->_enclosureHalfSize.resize(3);

    dump->_enclosureHalfSize[0] =  std::max(dump->_coreHalfSize[0] + dump->_minCoreShieldingThickness, enclosureHalfWidthMin);

    dump->_enclosureHalfSize[1] =  dump->_coreHalfSize[1] + dump->_minCoreShieldingThickness + dump->_magnetPitHalfSize[1];

    const double dumpPlusShieldingHalfLength =
      dump->_mouthHalfSize[2] + dump->_neutronCaveHalfSize[2] + dump->_coreHalfSize[2] + 0.5*dump->_minCoreShieldingThickness;

    const double totalFilterLength =
      dump->_collimator1._horizontalLength + 2*dump->_magnetPitHalfSize[2] + dump->_collimator2._horizontalLength;

    dump->_enclosureHalfSize[2] = std::max(dumpPlusShieldingHalfLength, 0.5*totalFilterLength);

    // shorthand notation
    const std::vector<double>& enclosureHalfSize = dump->_enclosureHalfSize;

    if(verbose) {
      std::cout<<"ProtonBeamDumpMaker"<<": ProtonBeamDump enclosure half size = ";
      std::copy(enclosureHalfSize.begin(), enclosureHalfSize.end(), std::ostream_iterator<double>(std::cout, ", "));
      std::cout<<std::endl;
    }

    // compute the position of the overall enclosure
    dump->_enclosureRotationInMu2e.rotateY(coreRotY);

    // The offset of the enclosure center w.r.t. the core, along the dump z
    const double coreOffset = 2*dump->_mouthHalfSize[2] + 2*dump->_neutronCaveHalfSize[2]
      + dump->_coreHalfSize[2] - enclosureHalfSize[2];

    dump->_enclosureCenterInMu2e[0] = dump->_coreCenterInMu2e[0] + coreOffset*sin(coreRotY);
    dump->_enclosureCenterInMu2e[1] = dump->_coreCenterInMu2e[1] + dump->_magnetPitHalfSize[1];
    dump->_enclosureCenterInMu2e[2] = dump->_coreCenterInMu2e[2] + coreOffset*cos(coreRotY);

    const CLHEP::Hep3Vector& enclosureCenterInMu2e = dump->_enclosureCenterInMu2e;

    // core relative to the enclosure
    dump->_coreCenterInEnclosure[0] = 0;
    dump->_coreCenterInEnclosure[1] = -dump->_magnetPitHalfSize[1];
    dump->_coreCenterInEnclosure[2] = -coreOffset;

    // position of the magnet pit
    dump->_magnetPitCenterInEnclosure[0] = 0.;
    dump->_magnetPitCenterInEnclosure[1] = enclosureHalfSize[1] - dump->_magnetPitHalfSize[1];
    dump->_magnetPitCenterInEnclosure[2] = enclosureHalfSize[2] - dump->_collimator1._horizontalLength - dump->_magnetPitHalfSize[2];

    // Shielding face coordinates
    dump->_shieldingFaceXmin = enclosureCenterInMu2e[0]
      + enclosureHalfSize[2] * sin(coreRotY)
      - enclosureHalfSize[0] * cos(coreRotY)
      ;

    dump->_shieldingFaceXmax = enclosureCenterInMu2e[0]
      + enclosureHalfSize[2] * sin(coreRotY)
      + enclosureHalfSize[0] * cos(coreRotY)
      ;

    dump->_shieldingFaceZatXmin = enclosureCenterInMu2e[2]
      + enclosureHalfSize[2] * cos(coreRotY)
      + enclosureHalfSize[0] * sin(coreRotY)
      ;

    dump->_shieldingFaceZatXmax = enclosureCenterInMu2e[2]
      + enclosureHalfSize[2] * cos(coreRotY)
      - enclosureHalfSize[0] * sin(coreRotY)
      ;

    //----------------------------------------------------------------
    // Compute the placement of filter elements

    // collimator1

    const CLHEP::Hep3Vector collimator1CenterInEnclosure(dump->coreCenterInEnclosure()[0]
                                                         + dump->filterEntranceOffsetX()
                                                         + 0.5*dump->collimator1().horizontalLength()*tan(dump->collimator1().angleH()),

                                                         dump->coreCenterInEnclosure()[1]
                                                         + dump->filterEntranceOffsetY()
                                                         + 0.5*dump->collimator1().horizontalLength()*tan(dump->collimator1().angleV())/cos(dump->filterAngleH()),

                                                         dump->enclosureHalfSize()[2] - 0.5*dump->collimator1().horizontalLength());

    dump->_collimator1CenterInEnclosure = collimator1CenterInEnclosure;

    {
      CLHEP::HepRotation tmp(CLHEP::HepRotation::IDENTITY);
      tmp.rotateX(entranceAngleV).rotateY(-angleH);
      dump->_collimator1RotationInMu2e = dump->_enclosureRotationInMu2e * tmp;
    }

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

    const double magnetAngleV =
      dump->filterEntranceAngleV() - dump->filterMagnet().trackBendHalfAngle(dump->extMonFilter_nominalMomentum());

    {
      CLHEP::HepRotation tmp(CLHEP::HepRotation::IDENTITY);
      tmp.rotateX(magnetAngleV).rotateY(-angleH);
      dump->_filterMagnetRotationInMu2e = dump->_enclosureRotationInMu2e * tmp;
    }

    // Arbitrary fix the (center,bottom,center) point of the magnet
    // aperture at the Z position of the magnet room center.

    // Z of magnet entrance at the (center,bottom) of aperture (in enclosure coords)
    const double magnetEntranceZ = dump->magnetPitCenterInEnclosure()[2]
      + dump->filterMagnet().outerHalfSize()[2] * cos(magnetAngleV) * cos(dump->filterAngleH());

    // Height of the entrance point (in enclosure coords)
    const double magnetEntranceY = dump->coreCenterInEnclosure()[1] + dump->filterEntranceOffsetY()
      + (dump->enclosureHalfSize()[2] - magnetEntranceZ)*tan(dump->filterEntranceAngleV())/cos(dump->filterAngleH())

      // the entrance point should be on the trajectory of reference
      // particle at the *bottom* of the collimator1 channel - shift for that:
      - 0.5 * dump->collimator1().channelHeight()[dump->collimator1().channelHeight().size()-1]/cos(dump->collimator1().angleV())
      ;

    // X of the entrance point (in enclosure coords)
    const double magnetEntranceX = dump->coreCenterInEnclosure()[0] + dump->filterEntranceOffsetX()
      + tan(dump->filterAngleH()) * (dump->enclosureHalfSize()[2] - magnetEntranceZ);

    // length from magnet entrance to magnet center in projection on (XZ)
    const double projLength =
      (dump->filterMagnet().outerHalfSize()[2]* cos(magnetAngleV) - 0.5*dump->filterMagnet().apertureHeight()*sin(magnetAngleV));

    // Compute the center of the magnet
    const CLHEP::Hep3Vector magnetCenterInEnclosure
      (magnetEntranceX + sin(dump->filterAngleH()) * projLength
       ,
       magnetEntranceY + sin(magnetAngleV)*dump->filterMagnet().outerHalfSize()[2] + cos(magnetAngleV)*0.5*dump->filterMagnet().apertureHeight()
       ,
       magnetEntranceZ - cos(dump->filterAngleH()) * projLength
       );

    dump->_filterMagnetCenterInEnclosure = magnetCenterInEnclosure;

    //----------------------------------------------------------------
    // collimator2: position so that the reference trajectory enters it at the bottom center
    // of the aperture

    const double col2CenterZ = -dump->_enclosureHalfSize[2] + 0.5*dump->_collimator2.horizontalLength();

    const double col2CenterX = dump->_coreCenterInEnclosure[0] + dump->_filterEntranceOffsetX +
      tan(dump->filterAngleH()) * (dump->_enclosureHalfSize[2] - col2CenterZ);

    // Z of the exit point from the magnet, at the (center,bottom) of aperture (in enclosure coords)
    const double magnetExitZ = dump->magnetPitCenterInEnclosure()[2]
      - dump->filterMagnet().outerHalfSize()[2] * cos(magnetAngleV) * cos(dump->filterAngleH());

    // Height of the exit point (in enclosure coords)
    const double magnetExitY = magnetEntranceY + sin(magnetAngleV) * 2*dump->_filterMagnet._outerHalfSize[2];

    const double col2CenterY = magnetExitY + (magnetExitZ - col2CenterZ)*tan(dump->_collimator2.angleV())/cos(dump->filterAngleH())
      + 0.5 * dump->_collimator2._channelHeight[0]/cos(dump->_collimator2.angleV());

    dump->_collimator2CenterInEnclosure = CLHEP::Hep3Vector(col2CenterX, col2CenterY, col2CenterZ);

    {
      CLHEP::HepRotation tmp(CLHEP::HepRotation::IDENTITY);
      tmp.rotateX(dump->_collimator2._angleV).rotateY(-angleH);
      dump->_collimator2RotationInMu2e = dump->_enclosureRotationInMu2e * tmp;
    }

    //----------------------------------------------------------------
    if(verbose) {
      std::cout<<"ProtonBeamDumpMaker"<<": Filter entrance offsets  = ("<<dump->_filterEntranceOffsetX<<", "<<dump->_filterEntranceOffsetY<<")"<<std::endl;
      std::cout<<"ProtonBeamDumpMaker"<<": coreOffset = "<<coreOffset<<std::endl;
      std::cout<<"ProtonBeamDumpMaker"<<": ProtonBeamDump core center in mu2e = "<<dump->_coreCenterInMu2e<<std::endl;
      std::cout<<"ProtonBeamDumpMaker"<<": ProtonBeamDump enclosure center in mu2e = "<<dump->_enclosureCenterInMu2e<<std::endl;
      std::cout<<"ProtonBeamDumpMaker"<<": coreCenterInEnclosure = "<<dump->_coreCenterInEnclosure<<std::endl;
      std::cout<<"ProtonBeamDumpMaker"<<": enclosure half size = { "
               <<dump->_enclosureHalfSize[0]<<", "
               <<dump->_enclosureHalfSize[1]<<", "
               <<dump->_enclosureHalfSize[2]<<" }"
               <<std::endl;
      std::cout<<"ProtonBeamDumpMaker"<<": magnetPitCenterInEnclosure = "<<dump->_magnetPitCenterInEnclosure<<std::endl;
      std::cout<<"ProtonBeamDumpMaker"<<": shieldingFaceXmin = "<<dump->_shieldingFaceXmin
               <<", Xmax = "<<dump->_shieldingFaceXmax<<std::endl;
      std::cout<<"ProtonBeamDumpMaker"<<": shieldingFaceZatXmin = "<<dump->_shieldingFaceZatXmin
               <<", ZatXmax = "<<dump->_shieldingFaceZatXmax<<std::endl;
      std::cout<<"ProtonBeamDumpMaker"<<": filter nominal momentum = "<<dump->extMonFilter_nominalMomentum()/CLHEP::GeV<<" GeV/c"<<std::endl;
      std::cout<<"ProtonBeamDumpMaker"<<": filterAngleH = "<<dump->filterAngleH()<<std::endl;
      std::cout<<"ProtonBeamDumpMaker"<<": filterAngleH in Mu2e, degrees= "<<(dump->coreRotY() - dump->filterAngleH())/CLHEP::degree<<std::endl;
      std::cout<<"ProtonBeamDumpMaker"<<": filter half bend angle  = "<<dump->filterMagnet().trackBendHalfAngle(dump->extMonFilter_nominalMomentum())<<std::endl;
      std::cout<<"ProtonBeamDumpMaker"<<": filter.angleV = "<<dump->filterEntranceAngleV()
               <<", c1.angleV  = "<<dump->collimator1().angleV()
               <<", magnet.angleV = "<<magnetAngleV
               <<", c2.angleV() = "<<dump->collimator2().angleV()<<std::endl;
      std::cout<<"ProtonBeamDumpMaker"<<": collimator1CenterInEnclosure = "<<dump->_collimator1CenterInEnclosure<<std::endl;
      std::cout<<"ProtonBeamDumpMaker"<<": collimator1.horizontalLength = "<<dump->_collimator1._horizontalLength<<std::endl;
      std::cout<<"ProtonBeamDumpMaker"<<": collimator2CenterInEnclosure = "<<dump->_collimator2CenterInEnclosure<<std::endl;
      std::cout<<"ProtonBeamDumpMaker"<<": collimator2.horizontalLength = "<<dump->_collimator2._horizontalLength<<std::endl;
    }

    if(verbose) {
      // Compute some positions in mu2e coordinates for verification and comparisons
      std::cout<<"ProtonBeamDumpMaker"<<": ProtonBeamDump::filterEntranceInMu2e() = "<<dump->filterEntranceInMu2e()<<std::endl;
      std::cout<<"ProtonBeamDumpMaker"<<" in Mu2e: col1 start = "
               <<enclosureToMu2ePoint(*dump.get(),
                                      // col1 start in enclosure coordinates
                                      dump->_collimator1CenterInEnclosure[0] - 0.5*dump->_collimator1._horizontalLength*tan(dump->_collimator1._angleH),
                                      dump->_collimator1CenterInEnclosure[1] - 0.5*dump->_collimator1._horizontalLength*tan(dump->_collimator1._angleV)/cos(dump->_filterAngleH),
                                      dump->_enclosureHalfSize[2]
                                      )
               <<std::endl;

      std::cout<<"ProtonBeamDumpMaker"<<" in Mu2e: col1 end = "
               <<enclosureToMu2ePoint(*dump.get(),
                                      dump->_collimator1CenterInEnclosure[0] + 0.5*dump->_collimator1._horizontalLength*tan(dump->_collimator1._angleH),
                                      dump->_collimator1CenterInEnclosure[1] + 0.5*dump->_collimator1._horizontalLength*tan(dump->_collimator1._angleV)/cos(dump->_filterAngleH),
                                      dump->_enclosureHalfSize[2] - dump->_collimator1._horizontalLength
                                      )
               <<std::endl;

      //----------------------------------------------------------------
      std::cout<<"ProtonBeamDumpMaker"<<" in Mu2e: col2 start = "
               <<enclosureToMu2ePoint(*dump.get(),
                                      dump->_collimator2CenterInEnclosure[0] - 0.5*dump->_collimator2._horizontalLength*tan(dump->_collimator2._angleH),
                                      dump->_collimator2CenterInEnclosure[1] - 0.5*dump->_collimator2._horizontalLength*tan(dump->_collimator2._angleV)/cos(dump->_filterAngleH),
                                      -dump->_enclosureHalfSize[2] + dump->_collimator2._horizontalLength
                                      )
               <<std::endl;

      std::cout<<"ProtonBeamDumpMaker"<<" in Mu2e: col2 end = "
               <<enclosureToMu2ePoint(*dump.get(),
                                      dump->_collimator2CenterInEnclosure[0] + 0.5*dump->_collimator2._horizontalLength*tan(dump->_collimator2._angleH),
                                      dump->_collimator2CenterInEnclosure[1] + 0.5*dump->_collimator2._horizontalLength*tan(dump->_collimator2._angleV)/cos(dump->_filterAngleH),
                                      -dump->_enclosureHalfSize[2]
                                      )
               <<std::endl;

      std::cout<<"ProtonBeamDumpMaker"<<": ProtonBeamDump::filterExitInMu2e() = "<<dump->filterExitInMu2e()<<std::endl;

      //----------------------------------------------------------------
      std::cout<<"ProtonBeamDumpMaker"<<": in Mu2e: magnet start = "
               <<enclosureToMu2ePoint(*dump.get(),
                                      dump->_filterMagnetCenterInEnclosure[0] - dump->_filterMagnet._outerHalfSize[2]*cos(magnetAngleV)*sin(dump->_filterAngleH),
                                      dump->_filterMagnetCenterInEnclosure[1] - dump->_filterMagnet._outerHalfSize[2]*sin(magnetAngleV),
                                      dump->_filterMagnetCenterInEnclosure[2] + dump->_filterMagnet._outerHalfSize[2]*cos(magnetAngleV)*cos(dump->_filterAngleH)
                                      )
               <<std::endl;

      std::cout<<"ProtonBeamDumpMaker"<<": in Mu2e: magnet end = "
               <<enclosureToMu2ePoint(*dump.get(),
                                      dump->_filterMagnetCenterInEnclosure[0] + dump->_filterMagnet._outerHalfSize[2]*cos(magnetAngleV)*sin(dump->_filterAngleH),
                                      dump->_filterMagnetCenterInEnclosure[1] + dump->_filterMagnet._outerHalfSize[2]*sin(magnetAngleV),
                                      dump->_filterMagnetCenterInEnclosure[2] - dump->_filterMagnet._outerHalfSize[2]*cos(magnetAngleV)*cos(dump->_filterAngleH)
                                      )
               <<std::endl;


    }

    return dump;

  } // make()

} // namespace mu2e
