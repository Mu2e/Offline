// Andrei Gaponenko, 2011

#include "Offline/GeometryService/inc/ExtMonFNALFilterMaker.hh"
#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALBuilding.hh"
#include "Offline/ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"

#include "Offline/GeometryService/inc/ExtMonFNALMagnetMaker.hh"

#include <algorithm>
#include <iterator>
#include <iostream>
#include <cmath>
#include <numeric>

#include "cetlib_except/exception.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include "Offline/ConfigTools/inc/SimpleConfig.hh"

//#define AGDEBUG(stuff) std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<": "<<stuff<<std::endl;
#define AGDEBUG(stuff)

namespace mu2e {

  namespace {
    std::vector<double> d2r(std::vector<double> v) {
      for(auto& x : v) { x/=2; }
      return v;
    }
  }

  //================================================================
  ExtMonFNALCollimator
  ExtMonFNALFilterMaker::readExtMonFNALCollimator(const std::string& name,
                                                  const SimpleConfig& c)
  {
    ExtMonFNALCollimator col;

    col._name = name;

    {
      std::vector<double> tmp;
      c.getVectorDouble("extMonFNAL."+name+".channelDiameter", tmp, 2);
      col._channelRadius = d2r(tmp);
    }

    {
      std::vector<double> tmp;
      c.getVectorDouble("extMonFNAL."+name+".alignmentPlugDiameter", tmp, 2);
      col._alignmentPlugRadius = d2r(tmp);
    }

    c.getVectorDouble("extMonFNAL."+name+".alignmentPlugInnerShellThickness", col._alignmentPlugInnerShellThickness, 2);
    c.getVectorDouble("extMonFNAL."+name+".alignmentPlugOuterShellThickness", col._alignmentPlugOuterShellThickness, 2);

    {
      std::vector<double> tmp;
      c.getVectorDouble("extMonFNAL."+name+".shotLinerInnerDiameter", tmp, 2);
      col._shotLinerInnerRadius = d2r(tmp);
    }

    c.getVectorDouble("extMonFNAL."+name+".shotLinerInnerThickness", col._shotLinerInnerThickness, 2);

    col._shotLinerOuterRadius = c.getDouble("extMonFNAL."+name+".shotLinerOuterDiameter")/2;
    col._shotLinerOuterThickness = c.getDouble("extMonFNAL."+name+".shotLinerOuterThickness");
    col._length = c.getDouble("extMonFNAL."+name+".length");
    col._radiusTransitiondZ = c.getDouble("extMonFNAL."+name+".radiusTransitiondZ");

    return col;
  }

  //================================================================
  ExtMonFNALFilter ExtMonFNALFilterMaker::read(const SimpleConfig& c,
                                               const ExtMonFNALBuilding& emfb,
                                               const ProtonBeamDump& dump)
  {
    using CLHEP::Hep3Vector;
    using CLHEP::Hep2Vector;

    ExtMonFNALFilter filter;

    const double pNominal = c.getDouble("extMonFNAL.filter.nominalMomentum") * CLHEP::MeV;
    filter.nominalMomentum_ = pNominal;

    //----------------------------------------------------------------
    const double angleH = c.getDouble("extMonFNAL.angleH") * CLHEP::radian;
    const double entranceAngleV = c.getDouble("extMonFNAL.entranceAngleV") * CLHEP::radian;

    filter.collimator1_ = readExtMonFNALCollimator("collimator1", c);
    filter.collimator1_.setFromDumpAngles(angleH, entranceAngleV, dump);

    //----------------------------------------------------------------
    // Position the entrance collimator

    const double filterEntranceOffsetX = c.getDouble("extMonFNAL.entranceOffsetX") * CLHEP::mm;
    const double filterEntranceOffsetY = c.getDouble("extMonFNAL.entranceOffsetY") * CLHEP::mm;

    // The offsets are for the point where the entrance collimator axis crosses
    // the "reference plane" parallel to the beam dump face (orthogonal to beam dump z)
    // at the distance dump.coreCenterDistanceToReferencePlane() from the core.

    const auto col1ReferenceInBeamDump = Hep3Vector(filterEntranceOffsetX,
                                                    filterEntranceOffsetY,
                                                    dump.coreCenterDistanceToReferencePlane());

    const auto col1ReferenceInMu2e = dump.beamDumpToMu2e_position(col1ReferenceInBeamDump);

    // The above posisitons the collimator axis in 2D, we also need to constrain
    // the position of the collimator along the axis.
    // The legacy convention is to have the collimator center 2.0 meters behind
    // the "reference plane".  (The length of the collimator used to be defined
    // by a 4.0 m thick concrete of the beam dump front.)

    const double referenceLength = dump.coreCenterDistanceToReferencePlane()
      - dump.coreCenterDistanceToShieldingFace()
      + dump.frontShieldingHalfSize()[2]; // FIXME: simplify - introduce a referenceLength parameter

    const auto col1DirInMu2e = filter.collimator1_.rotationInMu2e() * Hep3Vector(0,0,-1);
    const auto dumpDirInMu2e = dump.beamDumpToMu2e_momentum(Hep3Vector(0,0,-1));

    const double centerPar = referenceLength / dumpDirInMu2e.dot(col1DirInMu2e);
    const auto col1CenterInMu2e = col1ReferenceInMu2e + centerPar * col1DirInMu2e;

    filter.collimator1_._centerInMu2e = col1CenterInMu2e;

    //----------------------------------------------------------------
    // Filter magnet

    // the distance between the exit point of the reference trajectory from entrance collimator
    // and its entrance point to the magnet, on the magnet physical face.
    const double filterMagToColl = c.getDouble("extMonFNAL.filter.magnet.distanceToEntranceCollimator");

    const auto refTrajectoryEntranceInMu2e =
      filter.collimator1().trajectoryPointInMu2e(
                                                 -(filterMagToColl
                                                   + 0.5*filter.collimator1().length())
                                                 );

    filter.magnet_ = ExtMonFNALMagnetMaker::read(c,
                                                 "extMonFNAL.filter.magnet",
                                                 filter.collimator1().rotationInMu2e(),
                                                 refTrajectoryEntranceInMu2e,
                                                 pNominal);

    //----------------------------------------------------------------
    // Exit collimator

    filter.collimator2_ = readExtMonFNALCollimator("collimator2", c);

    // set the angles
    filter.collimator2_.setFromDumpAngles(angleH,
                                          entranceAngleV - 2 * filter.magnet().trackBendHalfAngle(pNominal),
                                          dump);

    // Set the position.  We put the center at the point where the
    // straight line of the nominal trajectory from the filter magnet
    // exit crosses the plane that is parallel to the collimator2
    // shielding wall and passes through the center of the wall.

    // unit vector along the nominal trajectory
    const CLHEP::Hep3Vector magDir = filter.magnet_.outRotationInMu2e()*CLHEP::Hep3Vector(0,0,-1);
    // point defining the nominal trajectory line
    const CLHEP::Hep3Vector magPoint = filter.magnet_.refPointInMu2e();

    // unit vector normal to the wall
    const CLHEP::Hep3Vector wallDir = dump.beamDumpToMu2e_momentum(CLHEP::Hep3Vector(0,0,-1));

    // Find a point with the Mu2e (x,z) coordinates in the middle of
    // the wall; the vertical (Mu2e Y) does not matter.

    // Find the 2D geometric center of the wall outline
    const auto coll2cog2d = std::accumulate(emfb.coll2ShieldingOutline().begin(),
                                            emfb.coll2ShieldingOutline().end(),
                                            Hep2Vector(0,0))
      / emfb.coll2ShieldingOutline().size();

    // Make it 3D by assigning Y=0, then shift to Mu2e
    const Hep3Vector wallPoint = Hep3Vector(coll2cog2d.y(), 0, coll2cog2d.x())
      + emfb.coll2ShieldingCenterInMu2e() ;

    filter.collimator2_._centerInMu2e  = magPoint + magDir*(wallDir.dot(wallPoint - magPoint)/wallDir.dot(magDir));

    //----------------------------------------------------------------
    int verbose = c.getInt("extMonFNAL.verbosityLevel");
    if(verbose) {
      std::cout<<"ExtMonFNALFilterMaker"<<": collimator1 angleH_inBeamDump= "<<filter.collimator1().angleH_inBeamDump()
               <<" ("  << filter.collimator1().angleH_inBeamDump() / CLHEP::degree <<" degree)"<<std::endl;

      std::cout<<"ExtMonFNALFilterMaker"<<": collimator1 angleH in Mu2e= "
               <<(dump.coreRotY() - filter.collimator1().angleH_inBeamDump())
               <<" ("  << (dump.coreRotY() - filter.collimator1().angleH_inBeamDump())/CLHEP::degree <<" degree)"<<std::endl;

      std::cout<<"ExtMonFNALFilterMaker"<<": collimator1 angleV= "<<filter.collimator1().angleV()
               <<" ("  << filter.collimator1().angleV() / CLHEP::degree <<" degree)"<<std::endl;

      //std::cout<<"ExtMonFNALFilterMaker"<<": collimator1 (dx/dz, dy/dz) = "<<filter.collimator1().dxdzdydz()<<std::endl;

      std::cout<<"ExtMonFNALFilterMaker"<<": collimator1 entranceInMu2e = "<<filter.collimator1().entranceInMu2e()<<std::endl;
      std::cout<<"ExtMonFNALFilterMaker"<<": collimator1 centerInMu2e = "<<filter.collimator1().centerInMu2e()<<std::endl;
      std::cout<<"ExtMonFNALFilterMaker"<<": collimator1 exitInMu2e = "<<filter.collimator1().exitInMu2e()<<std::endl;

      std::cout<<"ExtMonFNALFilterMaker"<<": ref traj entrace to filter magnet in Mu2e = "<< refTrajectoryEntranceInMu2e<<std::endl;
      std::cout<<"ExtMonFNALFilterMaker"<<": filter nominal momentum = "<<filter.nominalMomentum()/CLHEP::GeV<<" GeV/c"<<std::endl;
      std::cout<<"ExtMonFNALFilterMaker"<<": filter half bend angle  = "<<filter.magnet().trackBendHalfAngle(filter.magnet().nominalMomentum())<<std::endl;
      std::cout<<"ExtMonFNALFilterMaker"<<": filterMagnet().refPointInMu2e() = "<<filter.magnet().refPointInMu2e()<<std::endl;
      std::cout<<"ExtMonFNALFilterMaker"<<": filterMagnet().geometricCenterInMu2e() = "<<filter.magnet().geometricCenterInMu2e()<<std::endl;
      std::cout<<"ExtMonFNALFilterMaker"<<": filterMagnet().magnetRotationInMu2e() = "<<filter.magnet().magnetRotationInMu2e()<<std::endl;

      std::cout<<"ExtMonFNALFilterMaker"<<": filterEntranceInMu2e() = "<<filter.entranceInMu2e()<<std::endl;
      std::cout<<"ExtMonFNALFilterMaker"<<": filterExitInMu2e() = "<<filter.exitInMu2e()<<std::endl;
    }

    //----------------------------------------------------------------

    return filter;

  } // make()

} // namespace mu2e
