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

#define AGDEBUG(stuff) std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<": "<<stuff<<std::endl;
//#define AGDEBUG(stuff)

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

    //FIXME: tmp: ExtMonFNALFilter filter;
    ExtMonFNALFilter filter = emfb.filter(); // FIXME: use default ctr

    AGDEBUG("ORIG    col1 center = "<<filter.collimator1().centerInMu2e());
    AGDEBUG("ORIG    col1 rotation = "<<filter.collimator1().rotationInMu2e());

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

    AGDEBUG("col1ReferenceInMu2e = "<<col1ReferenceInMu2e);

    // The above posisitons the collimator axis in 2D, we also need to constrain
    // the position of the collimator along the axis.
    // The legacy convention is to have the collimator center 2.0 meters behind
    // the "reference plane".  (The length of the collimator used to be defined
    // by a 4.0 m thick concrete of the beam dump front.)

    const double referenceLength = dump.coreCenterDistanceToReferencePlane()
      - dump.coreCenterDistanceToShieldingFace()
      + dump.frontShieldingHalfSize()[2]; // FIXME: simplify - introduce a referenceLength parameter

    AGDEBUG("referenceLength = "<<referenceLength);

    const auto col1DirInMu2e = filter.collimator1_.rotationInMu2e() * Hep3Vector(0,0,-1);
    const auto dumpDirInMu2e = dump.beamDumpToMu2e_momentum(Hep3Vector(0,0,-1));

    AGDEBUG("col1DirInMu2e = "<<col1DirInMu2e);
    AGDEBUG("dumpDirInMu2e = "<<dumpDirInMu2e);

    const double centerPar = referenceLength / dumpDirInMu2e.dot(col1DirInMu2e);
    AGDEBUG("centerPar = "<<centerPar);

    const auto col1CenterInMu2e = col1ReferenceInMu2e + centerPar * col1DirInMu2e;

    AGDEBUG("NEW  col1 center = "<<col1CenterInMu2e);

    filter.collimator1_._centerInMu2e = col1CenterInMu2e;

//tmp:    //----------------------------------------------------------------
//tmp:    // Filter magnet
//tmp:
//tmp:    // the distance between the exit point of the reference trajectory from the upstream wall
//tmp:    // and its entrance to the magnet, on the magnet physical face.
//tmp:    const double filterMagToColl = c.getDouble("extMonFNAL.filter.magnet.distanceToUpstreamWall")
//tmp:      / (cos(angleH) * cos(entranceAngleV));
//tmp:
//tmp:    AGDEBUG("filterMagToColl computed = "<<filterMagToColl);
//tmp:    AGDEBUG("collimator1 exit in mu2e = "<<filter.collimator1().exitInMu2e());
//tmp:
//tmp:    const auto refTrajectoryEntranceInMu2e =
//tmp:      filter.collimator1().trajectoryPointInMu2e(
//tmp:                                                 -(filterMagToColl
//tmp:                                                   + 0.5*filter.collimator1().length())
//tmp:                                                 );
//tmp:
//tmp:    AGDEBUG("refTrajectoryEntranceInMu2e = "<<refTrajectoryEntranceInMu2e);
//tmp:
//tmp:    filter.magnet_ = ExtMonFNALMagnetMaker::read(c,
//tmp:                                                 "extMonFNAL.filter.magnet",
//tmp:                                                 filter.collimator1().rotationInMu2e(),
//tmp:                                                 refTrajectoryEntranceInMu2e,
//tmp:                                                 pNominal);
//tmp:
//tmp:    //----------------------------------------------------------------
//tmp:    // Exit collimator
//tmp:
//tmp:    filter.collimator2_ = readExtMonFNALCollimator("collimator2", c);
//tmp:
//tmp:    // set the angles
//tmp:    filter.collimator2_.setFromDumpAngles(angleH,
//tmp:                                          entranceAngleV - 2 * filter.magnet_.trackBendHalfAngle(pNominal),
//tmp:                                          dump);
//tmp:
//tmp:    // Set the position.  We put the center at the point where the
//tmp:    // straight line of the nominal trajectory from the filter magnet
//tmp:    // exit crosses the plane that is parallel to the collimator2
//tmp:    // shielding wall and passes through the center of the wall.
//tmp:
//tmp:    // unit vector along the nominal trajectory
//tmp:    const CLHEP::Hep3Vector magDir = filter.magnet_.outRotationInMu2e()*CLHEP::Hep3Vector(0,0,-1);
//tmp:    // point defining the nominal trajectory line
//tmp:    const CLHEP::Hep3Vector magPoint = filter.magnet_.refPointInMu2e();
//tmp:
//tmp:    // unit vector normal to the wall
//tmp:    const CLHEP::Hep3Vector wallDir = dump.beamDumpToMu2e_momentum(CLHEP::Hep3Vector(0,0,-1));
//tmp:
//tmp:    // Find a point with the Mu2e (x,z) coordinates in the middle of
//tmp:    // the wall; the vertical (Mu2e Y) does not matter.
//tmp:
//tmp:    // Find the 2D geometric center of the wall outline
//tmp:    const auto coll2cog2d = std::accumulate(emfb.coll2ShieldingOutline().begin(),
//tmp:                                            emfb.coll2ShieldingOutline().end(),
//tmp:                                            Hep2Vector(0,0))
//tmp:      / emfb.coll2ShieldingOutline().size();
//tmp:
//tmp:    // Make it 3D by assigning Y=0, then shift to Mu2e
//tmp:    const Hep3Vector wallPoint = Hep3Vector(coll2cog2d.y(), 0, coll2cog2d.x())
//tmp:      + emfb.coll2ShieldingCenterInMu2e() ;
//tmp:
//tmp:    filter.collimator2_._centerInMu2e  = magPoint + magDir*(wallDir.dot(wallPoint - magPoint)/wallDir.dot(magDir));
//tmp:
//tmp://    //----------------------------------------------------------------
//tmp://    int verbose = c.getInt("extMonFNAL.verbosityLevel");
//tmp://    if(verbose) {
//tmp://      std::cout<<"ExtMonFNALFilterMaker"<<": filter nominal momentum = "<<filter.nominalMomentum()/CLHEP::GeV<<" GeV/c"<<std::endl;
//tmp://      std::cout<<"ExtMonFNALFilterMaker"<<": filterAngleH = "<<emfb->filterAngleH()<<std::endl;
//tmp://      std::cout<<"ExtMonFNALFilterMaker"<<": filterAngleH in Mu2e, degrees= "<<(dump.coreRotY() - emfb->filterAngleH())/CLHEP::degree<<std::endl;
//tmp://      std::cout<<"ExtMonFNALFilterMaker"<<": filter half bend angle  = "<<emfb->filterMagnet().trackBendHalfAngle(emfb->filterMagnet().nominalMomentum())<<std::endl;
//tmp://      std::cout<<"ExtMonFNALFilterMaker"<<": filter.angleV = "<<emfb->filterEntranceAngleV()
//tmp://               <<", c1.angleV  = "<<emfb->collimator1().angleV()
//tmp://               <<", c2.angleV() = "<<emfb->collimator2().angleV()<<std::endl;
//tmp://      std::cout<<"ExtMonFNALFilterMaker"<<": collimator1CenterInMu2e = "<<emfb->_collimator1CenterInMu2e<<std::endl;
//tmp://      std::cout<<"ExtMonFNALFilterMaker"<<": collimator1 exit in Mu2e = "<< dump.beamDumpToMu2e_position(collimator1ExitInDump)<<std::endl;
//tmp://
//tmp://      std::cout<<"ExtMonFNALFilterMaker"<<": ref traj entrace to filter magnet in Mu2e = "<< dump.beamDumpToMu2e_position(refTrajFMEntranceInDump)<<std::endl;
//tmp://
//tmp://      std::cout<<"ExtMonFNALFilterMaker"<<": filterMagnet().refPointInMu2e() = "<<emfb->_filterMagnet.refPointInMu2e()<<std::endl;
//tmp://      std::cout<<"ExtMonFNALFilterMaker"<<": filterMagnet().geometricCenterInMu2e() = "<<emfb->_filterMagnet.geometricCenterInMu2e()<<std::endl;
//tmp://
//tmp://      std::cout<<"ExtMonFNALFilterMaker"<<": filterMagnet().magnetRotationInMu2e() = "<<emfb->_filterMagnet.magnetRotationInMu2e()<<std::endl;
//tmp://
//tmp://      std::cout<<"ExtMonFNALFilterMaker"<<": collimator2CenterInMu2e = "<<emfb->_collimator2CenterInMu2e<<std::endl;
//tmp://
//tmp://      std::cout<<"ExtMonFNALFilterMaker"<<": ExtMonFNALBuilding::filterEntranceInMu2e() = "<<emfb->filterEntranceInMu2e()<<std::endl;
//tmp://      std::cout<<"ExtMonFNALFilterMaker"<<": ExtMonFNALBuilding::filterExitInMu2e() = "<<emfb->filterExitInMu2e()<<std::endl;
//tmp://    }
//tmp:
//tmp:    //----------------------------------------------------------------

    return filter;

  } // make()

} // namespace mu2e
