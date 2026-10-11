//
// Construct and return Stopping Target Monitor (STM)
//
// Author: Anthony Palladino
// Update: Haichuan Cao Sept 2023
// Update: Yongyi Wu Sept-Oct 2026, hand-stacked downstream shielding
//         (stm.downstream.handstacked, STM_v11)
//
// Notes
// See mu2e-doc-XXXX for naming conventions etc.

// c++ includes
#include <iostream>
#include <iomanip>
#include <cmath>
#include <functional>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>

// clhep includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"

#include "Offline/GeneralUtilities/inc/OrientationResolver.hh"

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

// Mu2e includes
#include "Offline/GeometryService/inc/STMMaker.hh"
#include "Offline/STMGeom/inc/STM.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/Mu2eHallGeom/inc/Mu2eHall.hh"
#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/GeometryService/inc/G4GeometryOptions.hh"

using namespace std;

namespace mu2e {

  namespace {

    // Read a key that only a component being built needs. When it is
    // built the key is required as usual; when it is not, the geometry
    // file may leave the key out and it reads as zero or empty, since
    // nothing is placed from it.
    double readDoubleIfBuilt(SimpleConfig const & config,
                             std::string const & key, bool built) {
      return built ? config.getDouble(key) : config.getDouble(key, 0.);
    }
    bool readBoolIfBuilt(SimpleConfig const & config,
                         std::string const & key, bool built) {
      return built ? config.getBool(key) : config.getBool(key, false);
    }
    std::string readStringIfBuilt(SimpleConfig const & config,
                                  std::string const & key, bool built) {
      return built ? config.getString(key) : config.getString(key, "");
    }

    // How far a box reaches along one direction once it has been
    // turned.
    //
    // The hand-stacked sections describe a piece by its dimensions in
    // its own frame plus an orientation code, so those dimensions do
    // not say how much of any Mu2e axis the piece occupies until the
    // rotation is applied. Reaching for dx because a piece "looks
    // long in x" is the recurring mistake here: a 2x4x16 brick laid on
    // its side spans 2 in along a course whose dx says 16, and the
    // copper lining's 12.7 mm depth is its outline's v range under one
    // orientation and its sweep under another.
    //
    // dim is in the piece's own frame; rot is the placement rotation;
    // dir is the Mu2e direction asked about, and should be a unit
    // vector.
    //
    // constructSTM hands rot to G4PVPlacement as the frame rotation, so
    // the solid itself is turned by rot^-1: its own axis k lands on
    // rot^-1 * e_k in Mu2e. For a single 90 degree turn rot and rot^-1
    // put each axis on the same Mu2e line, so it makes no difference
    // there; for a two-axis code such as "101" they do not.
    double spanAlong(CLHEP::Hep3Vector const & dim,
                     CLHEP::HepRotation const & rot,
                     CLHEP::Hep3Vector const & dir) {
      double span = 0.;
      for (int k = 0; k < 3; ++k) {
        const CLHEP::Hep3Vector axis =
          rot.inverse() * CLHEP::Hep3Vector(k == 0, k == 1, k == 2);
        span += std::abs(axis.dot(dir)) * dim[k];
      }
      return span;
    }

    // The rotation an orientation code stands for.
    CLHEP::HepRotation rotationFor(std::string const & orientation) {
      CLHEP::HepRotation rot(CLHEP::HepRotation::IDENTITY);
      OrientationResolver OR;
      OR.getRotationFromOrientation(rot, orientation);
      return rot;
    }

    // How much of a brick type reaches along one Mu2e direction, with
    // the type range checked.
    //
    // dims is the list every section shares, indexed from 1 as the
    // geometry file numbers the types. `what` names the caller so a bad
    // type number says which section asked.
    double brickSpan(std::vector<CLHEP::Hep3Vector> const & dims,
                     int type,
                     CLHEP::HepRotation const & rot,
                     CLHEP::Hep3Vector const & dir,
                     std::string const & what) {
      if (type < 1 || type > int(dims.size())) {
        throw cet::exception("GEOM")
          << "STMMaker: " << what << " names brick type " << type
          << ", but stm.leadBrick.typeN defines only " << dims.size()
          << ".\n";
      }
      return spanAlong(dims[type-1], rot, dir);
    }

    // Read a section's bore list, stated as boreN then bore<i>Axis,
    // bore<i>R and bore<i>Offset<X/Y> under `base`, the section prefix
    // (e.g. "stm.FrontShieldingRight").
    //
    // The holes follow the collimator rather than the structure they
    // pass through, so each center comes from the spot offset with only
    // a correction of its own, and y is the beam plane.
    std::vector<BrickWallBore>
    readBores(SimpleConfig const & config, std::string const & base,
              double offsetSpot) {
      std::vector<BrickWallBore> bores;
      const int n = config.getInt(base + ".boreN");
      for (int i = 1; i <= n; ++i) {
        std::ostringstream key;
        key << base << ".bore" << i;
        BrickWallBore b;
        b.axis   = config.getString(key.str() + "Axis");
        b.radius = config.getDouble(key.str() + "R");
        b.offset = CLHEP::Hep3Vector(config.getDouble(key.str() + "OffsetX"),
                                     config.getDouble(key.str() + "OffsetY"), 0.);
        if (b.axis != "LaBr" && b.axis != "HPGe") {
          throw cet::exception("GEOM")
            << "STMMaker: " << key.str() << "Axis is \"" << b.axis
            << "\", which is not one of LaBr or HPGe.\n";
        }
        const double sign = (b.axis == "LaBr") ? +1. : -1.;
        b.center = CLHEP::Hep3Vector(sign*offsetSpot + b.offset.x(),
                                     b.offset.y(), 0.);
        bores.push_back(b);
      }
      return bores;
    }

    // The bore ids a piece names, range-checked against the section's
    // own bore list. `what` names the key so a bad id says which piece
    // asked for it.
    std::vector<int>
    readPieceBores(SimpleConfig const & config, std::string const & key,
                   size_t boreN, std::string const & what) {
      std::vector<int> ids;
      config.getVectorInt(key, ids, std::vector<int>());
      for (int id : ids) {
        if (id < 1 || id > int(boreN)) {
          throw cet::exception("GEOM")
            << "STMMaker: " << what << " names bore " << id
            << ", but " << boreN << " are defined for that section.\n";
        }
      }
      return ids;
    }

    // Which bores pass through the brick at a given position in a
    // course, and the wear to allow when checking they fit.
    //
    // A section that bores nothing leaves `bores` null and the whole
    // business is skipped. One that does bores passes its list and a
    // predicate saying which positions take which hole, and layCourse
    // checks each hole really lands on the brick named -- otherwise the
    // subtraction would clip nothing and the beam would see solid lead,
    // with nothing to say so.
    struct CourseBores {
      std::vector<BrickWallBore> const * bores;   // 0 when none
      // position in the course, 1-based -> bore ids through that brick
      std::function<std::vector<int>(int)>  at;
      double wearAlongCourse;
      double wearAlongPitch;
      std::string context;   // added to a failure, naming the caller

      CourseBores() : bores(0), wearAlongCourse(0.), wearAlongPitch(0.) {}
    };

    // Lay one course: butt a list of brick types end to end from
    // `start`, the corner the first brick's near face sits on, and hand
    // back the placed bricks. The three directions are BrickWall's.
    //
    // Each brick is centered half its own span inside the course along
    // all three, so the caller states a corner and gets centers back.
    // `step` returns the course's own span along the pitch and `thick`
    // its span along the depth: a course is one brick deep and one
    // brick tall.
    std::vector<BrickWallBrick>
    layCourse(std::vector<int> const & types,
              std::string const & orientation,
              CLHEP::Hep3Vector const & start,
              CLHEP::Hep3Vector const & courseDir,
              CLHEP::Hep3Vector const & pitchDir,
              CLHEP::Hep3Vector const & depthDir,
              std::vector<CLHEP::Hep3Vector> const & dims,
              std::string const & what,
              double * step = 0,
              double * thick = 0,
              CourseBores const & bores = CourseBores()) {
      if (types.empty()) {
        throw cet::exception("GEOM")
          << "STMMaker: " << what << " is an empty course.\n";
      }
      const CLHEP::HepRotation rot = rotationFor(orientation);

      const double h = brickSpan(dims, types[0], rot, pitchDir, what);
      const double d = brickSpan(dims, types[0], rot, depthDir, what);
      if (step)  *step  = h;
      if (thick) *thick = d;

      std::vector<BrickWallBrick> bricks;
      double along = 0.;
      for (size_t p = 0; p < types.size(); ++p) {
        const double w = brickSpan(dims, types[p], rot, courseDir, what);
        BrickWallBrick b;
        b.type        = types[p];
        b.orientation = orientation;
        b.center      = start
                      + courseDir*(along + w/2)
                      + pitchDir *(h/2)
                      + depthDir *(d/2);

        // The holes, if this section has any. They are pinned to the
        // beam while the course is pinned to the structure, so the two
        // can be moved apart: each is checked against the WORN brick,
        // which is the solid the hole is actually cut from.
        if (bores.bores) {
          for (int id : bores.at(int(p)+1)) {
            if (id < 1 || id > int(bores.bores->size())) {
              throw cet::exception("GEOM")
                << "STMMaker: " << bores.context << " names bore " << id
                << ", but only " << bores.bores->size() << " are defined.\n";
            }
            BrickWallBore const & bore = (*bores.bores)[id-1];

            const double hp = brickSpan(dims, types[p], rot, pitchDir, what);
            const double clo = b.center.dot(courseDir) - w/2 + bores.wearAlongCourse;
            const double chi = b.center.dot(courseDir) + w/2 - bores.wearAlongCourse;
            const double plo = b.center.dot(pitchDir)  - hp/2 + bores.wearAlongPitch;
            const double phi = b.center.dot(pitchDir)  + hp/2 - bores.wearAlongPitch;
            const double bc  = bore.center.dot(courseDir);
            const double bp  = bore.center.dot(pitchDir);

            if (bc - bore.radius < clo || bc + bore.radius > chi ||
                bp - bore.radius < plo || bp + bore.radius > phi) {
              throw cet::exception("GEOM")
                << "STMMaker: bore " << id << " (" << bore.axis
                << ") does not fit the brick at position " << p+1
                << " of " << what << ", which " << bores.context
                << " says is bored.\n"
                << "Along the course the hole spans ["
                << bc - bore.radius << ", " << bc + bore.radius
                << "] and the worn brick [" << clo << ", " << chi << "];\n"
                << "along the pitch the hole spans ["
                << bp - bore.radius << ", " << bp + bore.radius
                << "] and the worn brick [" << plo << ", " << phi << "] mm.\n"
                << "Check what pins the course against what pins the"
                << " beam, and this bore's own offset ("
                << bore.offset.x() << ", " << bore.offset.y() << " mm).\n";
            }
            b.bores.push_back(id);
          }
        }

        bricks.push_back(b);
        along += w;
      }
      return bricks;
    }

    // Lay one whole layer of lead: every course of it, butted side by
    // side along the pitch, and hand back the placed bricks.
    //
    // All four walls lay a layer the same way, so only the courses'
    // type lists are stated and where each lands follows from the
    // bricks. `prefix` is the layer's key stem, e.g.
    // "stm.BottomWall.leadLayer2", and <prefix>CourseN the count. Two
    // things vary between walls, and both fall back to a default:
    //
    //   <prefix>Course<c>FromRefZ     where the course starts along
    //                                 courseDir. Absent means 0 -- the
    //                                 side walls, whose courses all
    //                                 start at the layer corner.
    //   <prefix>Course<c>Orientation  read only when `orientation` is
    //                                 empty; the side walls pass one
    //                                 for the whole layer.
    //
    // `thick` returns the layer's depth, so the caller can advance its
    // cursor. A layer is one brick deep.
    std::vector<BrickWallBrick>
    layLayer(SimpleConfig const & config,
             std::string const & prefix,
             CLHEP::Hep3Vector const & start,
             CLHEP::Hep3Vector const & courseDir,
             CLHEP::Hep3Vector const & pitchDir,
             CLHEP::Hep3Vector const & depthDir,
             std::vector<CLHEP::Hep3Vector> const & dims,
             std::string const & orientation,
             double * thick = 0) {
      const int courseN = config.getInt(prefix + "CourseN");

      std::vector<BrickWallBrick> bricks;
      double acrossPitch = 0.;

      for (int c = 1; c <= courseN; ++c) {
        std::ostringstream ckey;
        ckey << prefix << "Course" << c;

        std::vector<int> types;
        config.getVectorInt(ckey.str(), types);

        const std::string orient = orientation.empty()
          ? config.getString(ckey.str() + "Orientation") : orientation;
        const double along = config.getDouble(ckey.str() + "FromRefZ", 0.);

        // The course's extent along the pitch comes back from layCourse
        // rather than being stated: it is what the next one starts
        // after.
        double extent = 0.;
        const std::vector<BrickWallBrick> course =
          layCourse(types, orient,
                    start + courseDir*along + pitchDir*acrossPitch,
                    courseDir, pitchDir, depthDir,
                    dims, ckey.str(), &extent, thick);
        bricks.insert(bricks.end(), course.begin(), course.end());
        acrossPitch += extent;
      }
      return bricks;
    }

    // The bricks of a layer that break its course pattern, each given
    // by its own center rather than by a place in a course:
    //
    //     {type, orientation, centerFromRefX, centerFromRefZ}
    //
    // The two offsets are Mu2e x and z as the geometry file states
    // them, NOT projections onto the layer's own directions -- a stray
    // is placed by measurement, so it is written the way it was
    // measured. Only the depth comes from the layer.
    //
    // The orientation rides in the tuple as a number, so it is turned
    // back into the three-digit code the resolver reads. Only the
    // bottom wall has strays, but the tuple form lives here so it is
    // written out once.
    std::vector<BrickWallBrick>
    layStrays(SimpleConfig const & config,
              std::string const & prefix,
              CLHEP::Hep3Vector const & start,
              CLHEP::Hep3Vector const & depthDir,
              std::vector<CLHEP::Hep3Vector> const & dims,
              double thick) {
      std::vector<BrickWallBrick> bricks;
      const int extraN = config.getInt(prefix + "ExtraN", 0);

      for (int e = 1; e <= extraN; ++e) {
        std::ostringstream ekey;
        ekey << prefix << "Extra" << e;
        std::vector<double> v;
        config.getVectorDouble(ekey.str(), v, 4);

        const int type = int(v[0]);
        if (type < 1 || type > int(dims.size())) {
          throw cet::exception("GEOM")
            << "STMMaker: " << ekey.str() << " names brick type " << type
            << ", but stm.leadBrick.typeN defines only " << dims.size()
            << ".\n";
        }
        // Leading zeros matter, hence the width.
        std::ostringstream orient;
        orient << std::setw(3) << std::setfill('0') << int(v[1]);

        BrickWallBrick b;
        b.type        = type;
        b.orientation = orient.str();
        b.center      = start + depthDir*(thick/2)
                      + CLHEP::Hep3Vector(v[2], 0., v[3]);
        bricks.push_back(b);
      }
      return bricks;
    }

  } // anonymous namespace

  // Constructor that gets information from the config file instead of
  // from arguments.
  STMMaker::STMMaker(SimpleConfig const & _config,
                     double solenoidOffset)
  {
    // if( ! _config.getBool("hasSTM",false) ) return;

    // create an empty STM
    _stm = unique_ptr<STM>(new STM());

    // access its object through a reference

    STM & stm = *_stm.get();

    parseConfig(_config);

    // One switch for the whole downstream geometry rather than a flag per
    // component; see STM::handstacked().
    stm._handstacked = _handstacked;

    // The shared standard lead bricks.
    if (_handstacked) {
      stm._pLeadBrickParams = std::unique_ptr<LeadBrick>
              (new LeadBrick(_leadBrickNames,
                             _leadBrickDims,
                             _leadBrickWearY,
                             _leadBrickWearXZ,
                             _leadBrickMaterial));
    }

    // now create the specific components

    // Fetch DS geom. object
    GeomHandle<DetectorSolenoid> ds;
    const CLHEP::Hep3Vector &dsP( ds->position() );

    //Create a reference position (most things in the STM geometry will be defined w.r.t. this position)
    //Our reference z is the downstream edge of the CRV-D as it was in crv_counters_v09.txt.
    //It is now set as a user variable in STM_v08.txt without changing its value.
    const CLHEP::Hep3Vector _STMMOffsetInMu2e(dsP.x(), 0.0, _stmReferenceZ );
    const CLHEP::HepRotation _magnetRotation = CLHEP::HepRotation::IDENTITY;
    double magnetZOffset = _magnetUpStrSpace+_magnetHalfLength;
    //calculate the magnet position assuming the shield pipe is flush to the wall
    if(_config.getBool("stm.magnet.usePipeAsOrigin", false)) {
      magnetZOffset = _magnetHalfLength + _shieldDnStrWallGap + 2.*_shieldPipeHalfLength + _shieldUpStrWallGap;
      if(!_shieldMatchPipeBlock) magnetZOffset += 2.*_shieldDnStrWallHalfLength;
    }
    const CLHEP::Hep3Vector _magnetOffsetInMu2e  = _STMMOffsetInMu2e + CLHEP::Hep3Vector(0.0,0.,magnetZOffset);
    const CLHEP::Hep3Vector _magnetHoleOffset  = CLHEP::Hep3Vector(_magnetHoleXOffset,_magnetHoleYOffset, 0.);

    //if (_magnetBuild){
      stm._pSTMMagnetParams = std::unique_ptr<PermanentMagnet>
        (new PermanentMagnet(_magnetBuild,
                             _magnetHalfWidth,
                             _magnetHalfHeight,
                             _magnetHalfLength,
                             _magnetHoleHalfWidth,
                             _magnetHoleHalfHeight,
                             _magnetOffsetInMu2e,
                             _magnetRotation,
                             _magnetHoleOffset,
                             _magnetMaterial,
                             _magnetHasLiner,
                             _magnetField,
                             _magnetFieldVisible
                            ));
    //}


    const CLHEP::HepRotation _FOVCollimatorRotation     = CLHEP::HepRotation::IDENTITY;
    CLHEP::Hep3Vector  _FOVCollimatorOffsetInMu2e = _magnetOffsetInMu2e + CLHEP::Hep3Vector(0.0,0.,_magnetHalfLength + _FOVCollimatorUpStrSpace + _FOVCollimatorHalfLength + _shieldPipeUpStrAirGap);
    // If we actually don't want to build the magnet, subtract off the offsets related to the magnet.
    // (We can't just set _magnetHalfLength = 0 in config because it is needed in various parts of constructSTM.cc
    // including to make a G4Box, which cannot have length 0...)
    if(!_magnetBuild) {
      _FOVCollimatorOffsetInMu2e -= CLHEP::Hep3Vector(0.0, 0.0, 2*_magnetHalfLength);
    }

    //if (_FOVCollimatorBuild){
      stm._pSTMFOVCollimatorParams = std::unique_ptr<STMCollimator>
        (new STMCollimator(_FOVCollimatorBuild,
                           _FOVCollimatorHalfWidth,
                           _FOVCollimatorHalfHeight,
                           _FOVCollimatorHalfLength,
                           _FOVCollimatorLinerBuild,
                           _FOVCollimatorLinerHalfWidth,
                           _FOVCollimatorLinerHalfHeight,
                           _FOVCollimatorLinerHalfLength,
                           _FOVCollimatorLinerCutOutHalfLength,
                           _FOVCollimatorHole1xOffset,
                           _FOVCollimatorHole1RadiusUpStr,
                           _FOVCollimatorHole1RadiusDnStr,
                           _FOVCollimatorHole1LinerBuild,
                           _FOVCollimatorHole1LinerThickness,
                           _FOVCollimatorHole2Build,
                           _FOVCollimatorHole2xOffset,
                           _FOVCollimatorHole2RadiusUpStr,
                           _FOVCollimatorHole2RadiusDnStr,
                           _FOVCollimatorHole2LinerBuild,
                           _FOVCollimatorHole2LinerThickness,
                           _FOVCollimatorOffsetInMu2e,
                           _FOVCollimatorRotation,
                           _FOVCollimatorMaterial,
                           _FOVCollimatorLinerMaterial,
                           _FOVCollimatorHoleLinerMaterial
                          ));
    //}


    const CLHEP::HepRotation _pipeRotation     = CLHEP::HepRotation::IDENTITY;
    const CLHEP::Hep3Vector  _pipeOffsetInMu2e = _STMMOffsetInMu2e + CLHEP::Hep3Vector(0.0,0.,_magnetUpStrSpace+_magnetHalfLength);
    //if (_pipeBuild){
      stm._pSTMTransportPipeParams = std::unique_ptr<TransportPipe>
        (new TransportPipe(_pipeBuild,
                           _pipeRadiusIn,
                           _pipeRadiusOut,
                           _pipeMaterial,
                           _pipeGasMaterial,
                           _pipeUpStrSpace,
                           _pipeDnStrHalfLength,
                           _pipeUpStrWindowMaterial,
                           _pipeUpStrWindowHalfLength,
                           _pipeDnStrWindowMaterial,
                           _pipeDnStrWindowHalfLength,
                           _pipeFlangeHalfLength,
                           _pipeFlangeOverhangR,
                           _pipeOffsetInMu2e,
                           _pipeRotation
                          ));
    //}

    double _magnetTableTopHalfWidth = 0.0;
    if ( _magnetBuild && !_FOVCollimatorBuild) _magnetTableTopHalfWidth = _magnetHalfWidth;
    if (!_magnetBuild &&  _FOVCollimatorBuild) _magnetTableTopHalfWidth = _FOVCollimatorHalfWidth;
    if ( _magnetBuild &&  _FOVCollimatorBuild) _magnetTableTopHalfWidth = std::max(_magnetHalfWidth,_FOVCollimatorHalfWidth);
    _magnetTableTopHalfWidth += _magnetTableTopExtraWidth;

    double _magnetTableTopHalfLength = 0.0;
    if ( _magnetBuild )        _magnetTableTopHalfLength += _magnetHalfLength;
    if ( _pipeBuild )          _magnetTableTopHalfLength += _pipeDnStrHalfLength;
    if ( _FOVCollimatorBuild ) _magnetTableTopHalfLength += 0.5*_FOVCollimatorUpStrSpace+_FOVCollimatorHalfLength;
    if ( _shieldBuild ) {
      _magnetTableTopHalfLength += 0.5*_shieldDnStrSpace+_shieldDnStrWallHalfLength;
      if (!_magnetBuild) {
        // Previous versions (<STM_v09) didn't need to include the shield pipe length in here.
        // Now we will include it, otherwise the table is tiny
        _magnetTableTopHalfLength += _shieldPipeHalfLength;
      }
    }
    _magnetTableTopHalfLength += _magnetTableTopExtraLength;

    const CLHEP::HepRotation _magnetTableRotation     = CLHEP::HepRotation::IDENTITY;
    CLHEP::Hep3Vector  _magnetTableOffsetInMu2e = _STMMOffsetInMu2e - CLHEP::Hep3Vector(0.0,_magnetHalfHeight+_magnetTableTopHalfHeight,0.0);
    if ( _magnetBuild )        _magnetTableOffsetInMu2e += CLHEP::Hep3Vector(0.0,0.,_magnetUpStrSpace+_magnetHalfLength);
    if ( _pipeBuild )          _magnetTableOffsetInMu2e += CLHEP::Hep3Vector(0.0,0.,_pipeDnStrHalfLength);
    if ( _FOVCollimatorBuild ) _magnetTableOffsetInMu2e += CLHEP::Hep3Vector(0.0,0.,0.5*_FOVCollimatorUpStrSpace+_FOVCollimatorHalfLength);
    if ( _shieldBuild ) {
      _magnetTableOffsetInMu2e += CLHEP::Hep3Vector(0.0,0.,-0.5*_shieldDnStrSpace-_shieldDnStrWallHalfLength + _shieldPipeUpStrAirGap);
      if (!_magnetBuild) {
        // Previous versions (<STM_v09) didn't need to include the shield pipe length in here.
        // Now we will include it, otherwise the table is tiny
        _magnetTableOffsetInMu2e += CLHEP::Hep3Vector(0.0,0.,_shieldPipeHalfLength+0.5*_shieldDnStrSpace+_shieldDnStrWallHalfLength+_shieldPipeUpStrAirGap);
      }
    }


    //if (_magnetTableBuild && (_magnetBuild||_FOVCollimatorBuild) ){
      stm._pSTMMagnetSupportTableParams = std::unique_ptr<SupportTable>
        (new SupportTable( _magnetTableBuild,
                           _magnetTableTopHalfWidth,
                           _magnetTableTopHalfHeight,
                           _magnetTableTopHalfLength,
                           _magnetTableLegRadius,
                           _magnetTableOffsetInMu2e,
                           _magnetTableRotation,
                           _magnetTableMaterial
                          ));
    //}


      ////////////////////////////////////
      // STM Downstream Area
      //
    //The STM geometry must fit inside the detector hall, so find the z of the East hall wall
    GeomHandle<Mu2eHall> hall;
    const double z_hall_inside_max = hall->getWallExtentz("dsArea",1)/CLHEP::mm;//the integer allows you to specify which side of which wall you want the z for: 1 = west side of east wall (i.e. the z of the inside surface of the east wall)
    const CLHEP::Hep3Vector BeamAxisAtEastWallInMu2e(dsP.x(), 0.0, z_hall_inside_max );
    const double yExtentLow = std::abs(_config.getDouble("yOfFloorSurface.below.mu2eOrigin") );
    const CLHEP::Hep3Vector FloorAtEastWallInMu2e = BeamAxisAtEastWallInMu2e - CLHEP::Hep3Vector(0.0, yExtentLow, 0.0);

    // Define the envelope w.r.t the floor at the east wall
    const CLHEP::HepRotation _stmDnStrEnvRotation     = CLHEP::HepRotation::IDENTITY;
    const CLHEP::Hep3Vector _stmDnStrEnvPositionInMu2e = FloorAtEastWallInMu2e + CLHEP::Hep3Vector(0.0, +_stmDnStrEnvHalfHeight, -_stmDnStrEnvHalfLength);
    stm._pSTMDnStrEnvParams = std::unique_ptr<STMDownstreamEnvelope>
      (new STMDownstreamEnvelope(_stmDnStrEnvBuild,
                                 _stmDnStrEnvHalfWidth,
                                 _stmDnStrEnvHalfHeight,
                                 _stmDnStrEnvHalfLength,
                                 _stmDnStrEnvPositionInMu2e,
                                 _stmDnStrEnvRotation,
                                 _stmDnStrEnvMaterial
                                 ));

    const CLHEP::HepRotation _SSCollimatorRotation     = CLHEP::HepRotation::IDENTITY;
    const CLHEP::Hep3Vector  _SSCollimatorOffsetInMu2e = BeamAxisAtEastWallInMu2e + CLHEP::Hep3Vector(0.0,0.,-_stmZAllowed+_SSCollimatorHalfLength);
    //if(_SSCollimatorBuild){
      stm._pSTMSSCollimatorParams = std::unique_ptr<STMCollimator>
        (new STMCollimator(_SSCollimatorBuild,
                           _SSCollimatorHalfWidth,
                           _SSCollimatorHalfHeight,
                           _SSCollimatorHalfLength,
                           _SSCollimatorLinerBuild,
                           _SSCollimatorLinerHalfWidth,
                           _SSCollimatorLinerHalfHeight,
                           _SSCollimatorLinerHalfLength,
                           _SSCollimatorLinerCutOutHalfLength,
                           _SSCollimatorHole1xOffset,
                           _SSCollimatorHole1RadiusUpStr,
                           _SSCollimatorHole1RadiusDnStr,
                           _SSCollimatorHole1LinerBuild,
                           _SSCollimatorHole1LinerThickness,
                           _SSCollimatorHole2Build,
                           _SSCollimatorHole2xOffset,
                           _SSCollimatorHole2RadiusUpStr,
                           _SSCollimatorHole2RadiusDnStr,
                           _SSCollimatorHole2LinerBuild,
                           _SSCollimatorHole2LinerThickness,
                           _SSCollimatorOffsetInMu2e,
                           _SSCollimatorRotation,
                           _SSCollimatorMaterial,
                           _SSCollimatorLinerMaterial,
                           _SSCollimatorHoleLinerMaterial
                          ));
    //}


    double _detectorTableTopHalfWidth = _SSCollimatorHalfWidth + _detectorTableTopExtraWidth;
    double _detectorTableTopHalfLength = 0.5*_stmZAllowed - 1.0;
    const CLHEP::HepRotation _detectorTableRotation = CLHEP::HepRotation::IDENTITY;
    CLHEP::Hep3Vector  _detectorTableOffsetInMu2e = BeamAxisAtEastWallInMu2e + CLHEP::Hep3Vector(0.0,-_SSCollimatorHalfHeight-_detectorTableTopHalfHeight, -_stmZAllowed+_detectorTableTopHalfLength);

    //if (_detectorTableBuild && _SSCollimatorBuild ){
      stm._pSTMDetectorSupportTableParams = std::unique_ptr<SupportTable>
        (new SupportTable( _detectorTableBuild,
                           _detectorTableTopHalfWidth,
                           _detectorTableTopHalfHeight,
                           _detectorTableTopHalfLength,
                           _detectorTableLegRadius,
                           _detectorTableOffsetInMu2e,
                           _detectorTableRotation,
                           _detectorTableMaterial
                          ));
    //}


    const CLHEP::HepRotation _detector1Rotation = CLHEP::HepRotation::IDENTITY;
    const CLHEP::Hep3Vector  _detector1OffsetInMu2e  = _SSCollimatorOffsetInMu2e + CLHEP::Hep3Vector(_detector1xOffset, 0.0, _SSCollimatorHalfLength+_detector1CanUpStrSpace+_detector1CanHalfLength);
    //if (_detector1Build){
      stm._pSTMDetector1Params = std::unique_ptr<GeDetector>
        (new GeDetector(_detector1Build,
                        _detector1CrystalMaterial,
                        _detector1CrystalRadiusIn,
                        _detector1CrystalRadiusOut,
                        _detector1CrystalHalfLength,
                        _detector1CanMaterial,
                        _detector1CanRadiusIn,
                        _detector1CanRadiusOut,
                        _detector1CanHalfLength,
                        _detector1CanUpStrWindowMaterial,
                        _detector1CanUpStrWindowHalfLength,
                        _detector1CanGasMaterial,
                        _detector1OffsetInMu2e,
                        _detector1Rotation
                       ));
      //}

    const CLHEP::HepRotation _detector2Rotation = CLHEP::HepRotation::IDENTITY;
    const CLHEP::Hep3Vector  _detector2OffsetInMu2e  = _SSCollimatorOffsetInMu2e + CLHEP::Hep3Vector(_detector2xOffset, 0.0, _SSCollimatorHalfLength+_detector2CanUpStrSpace+_detector2CanHalfLength);
    //if (_detector2Build){
      stm._pSTMDetector2Params = std::unique_ptr<GeDetector>
        (new GeDetector(_detector2Build,
                        _detector2CrystalMaterial,
                        _detector2CrystalRadiusIn,
                        _detector2CrystalRadiusOut,
                        _detector2CrystalHalfLength,
                        _detector2CanMaterial,
                        _detector2CanRadiusIn,
                        _detector2CanRadiusOut,
                        _detector2CanHalfLength,
                        _detector2CanUpStrWindowMaterial,
                        _detector2CanUpStrWindowHalfLength,
                        _detector2CanGasMaterial,
                        _detector2OffsetInMu2e,
                        _detector2Rotation
                       ));
      //}


    const CLHEP::HepRotation _shieldRotation = CLHEP::HepRotation::IDENTITY;
    const CLHEP::Hep3Vector  _shieldOffsetInMu2e  = _FOVCollimatorOffsetInMu2e + CLHEP::Hep3Vector(0.0, 0.0, -_FOVCollimatorHalfLength);
    //if (_shieldBuild){
      stm._pSTMShieldPipeParams = std::unique_ptr<ShieldPipe>
        (new ShieldPipe(_shieldBuild,
                        _shieldRadiusIn,
                        _shieldHasLiner,
                        _shieldLinerWidth,
                        _shieldRadiusOut,
                        _shieldPipeHalfLength,
                        _shieldMaterialLiner,
                        _shieldMaterial,
                        _shieldMatchPipeBlock,
                        _shieldUpStrSpace,
                        _shieldDnStrSpace,
                        _shieldDnStrWallHalfLength,
                        _shieldDnStrWallHoleRadius,
                        _shieldDnStrWallHalfHeight,
                        _shieldDnStrWallHalfWidth,
                        _shieldDnStrWallGap,
                        _shieldDnStrWallMaterial,
                        _shieldBuildMatingBlock,
                        _shieldPipeUpStrAirGap,
                        _shieldOffsetInMu2e, //This is upstream edge of FOV collimator for now.
                        _shieldRotation
                       ));
    //}

   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   /// The geometries below were updated by Haichuan Cao in Sept. 2023
   /// and by Yongyi Wu in Sept-Oct 2026 for the hand-stacked downstream
   /// shielding (stm.downstream.handstacked).

    // NOT the point the downstream STM is placed from. constructSTM places
    // the SSC, its support, the shield house, the detectors and the back
    // and electronic shielding from its own STMShieldingRef, which is
    // BackShielding BPThick + ShieldingPipeGap further upstream than this
    // (east wall inner face - FrontToWall - BPThick - ShieldingPipeGap).
    //
    // Of the origins built from this point below, only STM_SSC's is read,
    // and only as a carrier: constructSTM and VirtualDetectorMaker take it,
    // subtract Wdepth_f/2 and BPThick + ShieldingPipeGap, and so recover
    // the reference above. The others (SSCSupport, FrontShielding,
    // HPGeDetector, LaBrDetector, Bottom/Left/Right/TopShielding,
    // BackShielding) are stored but not read, and do not match where those
    // pieces are placed.
    const CLHEP::Hep3Vector  _STMShieldingRef = BeamAxisAtEastWallInMu2e + CLHEP::Hep3Vector(0., 0., -_STM_SSCFrontToWall);

   ////////////////////////////////////////////////////////////////
   //Tungsten Spot-Size Collimator

    const CLHEP::HepRotation _STM_SSCRotation     = CLHEP::HepRotation::IDENTITY;
    const CLHEP::Hep3Vector  _STM_SSCOffsetInMu2e = _STMShieldingRef + CLHEP::Hep3Vector(0, 0, _STM_SSCWdepth_f/2);

    if (_handstacked) {
      // The updated SSC: a single symmetric block, bores stepped in z.
          stm._pSTM_SSCParams = std::unique_ptr<STM_SSC>
                  (new STM_SSC(_STM_SSCBuild,
                               _STM_SSCVDBuild,
                               _STM_SSCW_width,
                               _STM_SSCW_height,
                               _STM_SSCWdepth_f,
                               _STM_SSCWdepth_b,
                               _STM_SSCr_LaBr_f,
                               _STM_SSCr_HPGe_f,
                               _STM_SSCr_LaBr_b,
                               _STM_SSCr_HPGe_b,
                               _STM_SSCoffset_Spot,
                               _STM_SSCleak,
                               _STM_SSCFrontToWall,
                               _STM_SSCZGap,
                               _STM_SSCZGapBack,
                               _STM_SSCboreToBase,
                               _STM_SSCOffsetInMu2e,
                               _STM_SSCRotation,
                               _STM_SSCMaterial));
    } else {
          stm._pSTM_SSCParams = std::unique_ptr<STM_SSC>
                  (new STM_SSC(_STM_SSCBuild,
                               _STM_SSCVDBuild,
                               _STM_SSCdelta_WlR,
                               _STM_SSCdelta_WlL,
                               _STM_SSCW_middle,
                               _STM_SSCW_height,
                               _STM_SSCWdepth_f,
                               _STM_SSCWdepth_b,
                               _STM_SSCAperture_HPGe1,
                               _STM_SSCAperture_HPGe2,
                               _STM_SSCAperture_LaBr1,
                               _STM_SSCAperture_LaBr2,
                               _STM_SSCoffset_Spot,
                               _STM_SSCleak,
                               _STM_SSCFrontToWall,
                               _STM_SSCZGap,
                               _STM_SSCZGapBack,
                               _STM_SSCOffsetInMu2e,
                               _STM_SSCRotation,
                               _STM_SSCMaterial));
    }

   ////////////////////////////////////////////////////////////////
   //Spot-Size Collimator Support

    const CLHEP::HepRotation _SSCSupportRotation     = CLHEP::HepRotation::IDENTITY;
    const CLHEP::Hep3Vector  _SSCSupportOffsetInMu2e = _STMShieldingRef + CLHEP::Hep3Vector(0, 0, _STM_SSCZGap/2);

    if (_handstacked) {
      // The updated support: a five-plate steel cradle around the collimator.
          stm._pSSCSupportParams = std::unique_ptr<SSCSupport>
                  (new SSCSupport(_SSCSupportBuild,
                               _SSCSupportdepth,
                               _SSCSupportside_T,
                               _SSCSupportside_H,
                               _SSCSupportplate_base_T,
                               _SSCSupportbottom_T,
                               _SSCSupporttop_T,
                               _SSCSupportMaterial,
                               _SSCSupportOffsetInMu2e,
                               _SSCSupportRotation));
    } else {
          stm._pSSCSupportParams = std::unique_ptr<SSCSupport>
                  (new SSCSupport(_SSCSupportBuild,
                               _SSCSupporttable_L,
                               _SSCSupporttable_H,
                               _SSCSupporttable_T,
                               _SSCSupportleg_L,
                               _SSCSupportleg_H,
                               _SSCSupportleg_T,
                               _SSCSupportbase_L,
                               _SSCSupportbase_H,
                               _SSCSupportbase_T,
                               _SSCSupportwall_L,
                               _SSCSupportwall_H,
                               _SSCSupportwall_T,
                               _SSCSupporthole_H,
                               _SSCSupporthole_T,
                               _SSCSupportFLeadStand_L,
                               _SSCSupportFLeadStand_H,
                               _SSCSupportFLeadStand_T,
                               _SSCSupportFLeadShim_H,
                               _SSCSupportFLeadShim_T,
                               _SSCSupportFAluminumShim_T,
                               _SSCSupportFAluminumExtra_L,
                               _SSCSupportFAluminumExtra_H,
                               _SSCSupportOffsetInMu2e,
                               _SSCSupportRotation));
    }

   ////////////////////////////////////////////////////////////////
   //STM Front Shielding
   //
   // The updated geometry splits this into a left and a right half.
   // The right one -- the section behind the collimator -- is built
   // here in place of the earlier single description; its lead courses
   // and poly sheets were expanded from the config's type lists in
   // parseConfig, so what is passed on are finished placements.

    // The earlier front shielding's depth and width, used further down
    // only by the earlier Bottom, Left, Right and Top sections. The
    // updated walls replace those and place themselves off the right
    // half instead, so these stay zero under the updated geometry and are
    // set only by the branch that needs them.
    double _FrontS_Thickness = 0.;
    double _FrontS_Length    = 0.;

    if (_handstacked) {
      stm._pFrontShieldingRightParams = std::unique_ptr<FrontShieldingRight>
              (new FrontShieldingRight(_FrontShieldingRightBuild,
                                       _FrontShieldingRightLeadLayers,
                                       _FrontShieldingRightSheets,
                                       _FrontShieldingRightBlockMaterial,
                                       _FrontShieldingRightBlockHalfDim,
                                       _FrontShieldingRightBlockCenter,
                                       _FrontShieldingRightPipes,
                                       _FrontShieldingRightPlate,
                                       _FrontShieldingRightBackZ));

      stm._pFrontShieldingLeftParams = std::unique_ptr<FrontShieldingLeft>
              (new FrontShieldingLeft(_FrontShieldingLeftBuild,
                                      _FrontShieldingLeftSheets,
                                      _FrontShieldingLeftBrickGroups,
                                      _FrontShieldingLeftPrism));
    } else {

      const CLHEP::HepRotation _FrontSRotation   = CLHEP::HepRotation::IDENTITY;
      _FrontS_Thickness = _FrontStungstendepth + _FrontSLeakForSSC + _FrontSleaddepth2*3 + _FrontSBPdepth*2 + _FrontScopperdepth;
      _FrontS_Length    = _FrontStungstenlength + 2*_FrontSLeakForSSC + _FrontSfPb_lengthL + _FrontSfPb_lengthR;

      const CLHEP::Hep3Vector  _FrontSOffsetInMu2e = _STMShieldingRef + CLHEP::Hep3Vector(0, 0, _FrontS_Thickness/2);

          stm._pSTMFrontShieldingParams = std::unique_ptr<FrontShielding>
          (new FrontShielding(_FrontShieldingBuild,
                              _FrontSHeightofRoom,
                              _FrontStungstenlength,
                              _FrontStungstendepth,
                              _FrontSleaddepth1,
                              _FrontSleaddepth2,
                              _FrontSaluminumdepth,
                              _FrontScopperdepth,
                              _FrontSBPdepth,
                              _FrontSfPb_lengthL,
                              _FrontSfPb_lengthR,
                              _FrontSGapForTop,
                              _FrontSLeakForSSC,
                              _FrontSCopperL,
                              _FrontS_H,
                              _FrontSHole_r,
                              _FrontSOffsetInMu2e,
                              _FrontSRotation));
    }

   ////////////////////////////////////////////////////////////////
   //STM HPGe Detector

    CLHEP::HepRotation rotHPGe = CLHEP::HepRotation::IDENTITY;
    rotHPGe.rotateY(45*CLHEP::degree);

    const CLHEP::HepRotation _HPGeRotation     =  rotHPGe;
    const CLHEP::Hep3Vector  _HPGeOffsetInMu2e = _STMShieldingRef + CLHEP::Hep3Vector(-_STM_SSCoffset_Spot + _HPGeoffset_HPGe, 0., _HPGeZ_HPGe);
          stm._pSTMHPGeDetectorParams = std::unique_ptr<HPGeDetector>
          (new HPGeDetector(_HPGeBuild,
                            _HPGecrystalMaterial,
                            _HPGeholeMaterial,
                            _HPGewindowMaterial,
                            _HPGewallMaterial,
                            _HPGecapsuleMaterial,
                            _HPGeEndcapR,
                            _HPGeEndcapL,
                            _HPGeCrystalR,
                            _HPGeCrystalL,
                            _HPGeZ_HPGe,
                            _HPGeHoleR,
                            _HPGeHoleL,
                            _HPGeCapsule_Wallthick,
                            _HPGeCapsule_Windowthick,
                            _HPGeCapsule_Endthick,
                            _HPGeCapsule_Walllength,
                            _HPGeWindowD,
                            _HPGeEndcapD,
                            _HPGeAirD,
                            _HPGeoffset_HPGe,
                            _HPGeOffsetInMu2e,
                            _HPGeRotation));


   ////////////////////////////////////////////////////////////////
   //STM LaBr Detector

    const CLHEP::HepRotation _LaBrRotation     = CLHEP::HepRotation::IDENTITY;
    const CLHEP::Hep3Vector  _LaBrOffsetInMu2e = _STMShieldingRef + CLHEP::Hep3Vector(_STM_SSCoffset_Spot + _LaBroffset_LaBr, 0., _LaBrZ_LaBr);
          stm._pSTMLaBrDetectorParams = std::unique_ptr<LaBrDetector>
          (new LaBrDetector(_LaBrBuild,
                            _LaBrcrystalMaterial,
                            _LaBrwindowMaterial,
                            _LaBrwallMaterial,
                            _LaBrEndcapR,
                            _LaBrEndcapL,
                            _LaBrCrystalR,
                            _LaBrCrystalL,
                            _LaBrZ_LaBr,
                            _LaBrWindowD,
                            _LaBrEndcapD,
                            _LaBrAirD,
                            _LaBroffset_LaBr,
                            _LaBrOffsetInMu2e,
                            _LaBrRotation));


   ////////////////////////////////////////////////////////////////
   //STM Bottom Shielding

    if (_handstacked) {
      // The updated floor: the baseplate plus five layers worked upward
      // from it, already placed in parseConfig.
          stm._pBottomWallParams = std::unique_ptr<BottomWall>
          (new BottomWall(_BottomWallBuild,
                          _BottomWallPlate,
                          _BottomWallLeadLayers,
                          _BottomWallPrisms));
    } else {

      const double _BottomS_Thickness = _BottomSleaddepth*2 + _BottomSBPdepth*2 + _BottomScopperdepth;

      const double B_dX = - _FrontS_Length - _BottomSfloor_Zlength/4 + 264.12;
      const double B_dY = -_FrontSHeightofRoom/2 - _BottomS_Thickness/2;
      const double B_dZ = _FrontS_Thickness - 12.7 + _BottomSfloor_Zlength/2;

      const CLHEP::HepRotation _BottomSRotation   = CLHEP::HepRotation::IDENTITY;
      const CLHEP::Hep3Vector  _BottomSOffsetInMu2e = _STMShieldingRef + CLHEP::Hep3Vector(B_dX, B_dY, B_dZ);


          stm._pSTMBottomShieldingParams = std::unique_ptr<BottomShielding>
          (new BottomShielding(_BottomShieldingBuild,
                               _BottomSfloor_Zlength,
                               _BottomSFront_LB,
                               _BottomSFront_LB_inner,
                               _BottomSleaddepth,
                               _BottomScopperdepth,
                               _BottomSBPdepth,
                               _BottomSOffsetInMu2e,
                               _BottomSRotation));
    }


   ////////////////////////////////////////////////////////////////
   //STM Left Shielding

    if (_handstacked) {
      // The updated left wall: five layers worked inward from the
      // outermost, already placed in parseConfig.
          stm._pLeftWallParams = std::unique_ptr<LeftWall>
          (new LeftWall(_LeftWallBuild,
                        _LeftWallLeadLayers,
                        _LeftWallSheets));
    } else {

      const double _LeftS_Thickness = _LeftSleaddepth*2 + _LeftSBPdepth*2 + _LeftScopperdepth;

      const double L_dX = _LeftSXmin + _LeftS_Thickness/2;
      const double L_dY = 0;
      const double L_dZ = _FrontS_Thickness + _LeftS_Length/2;

      const CLHEP::HepRotation _LeftSRotation   = CLHEP::HepRotation::IDENTITY;
      const CLHEP::Hep3Vector  _LeftSOffsetInMu2e = _STMShieldingRef + CLHEP::Hep3Vector(L_dX, L_dY, L_dZ);

          stm._pSTMLeftShieldingParams = std::unique_ptr<LeftShielding>
          (new LeftShielding(_LeftShieldingBuild,
                             _LeftS_Length,
                             _LeftSleaddepth,
                             _LeftScopperdepth,
                             _LeftSBPdepth,
                             _LeftSXmin,
                             _LeftSOffsetInMu2e,
                             _LeftSRotation));
    }


   ////////////////////////////////////////////////////////////////
   //STM Right Shielding

    if (_handstacked) {
      // The updated right wall: five layers worked inward from the
      // front shielding, already placed in parseConfig.
          stm._pRightWallParams = std::unique_ptr<RightWall>
          (new RightWall(_RightWallBuild,
                         _RightWallLeadLayers,
                         _RightWallSheets,
                         _RightWallEdgePrism));
    } else {

      const double _RightS_Thickness = (_RightSleaddepth*2 + _RightSBPdepth*2 + _RightScopperdepth)*sqrt(2);

      const double R_dX = _RightSXmax - _RightS_Thickness/2;
      const double R_dY = 0;
      const double R_dZ = _FrontS_Thickness + _RightS_Length/2;

      const CLHEP::HepRotation _RightSRotation   = CLHEP::HepRotation::IDENTITY;
      const CLHEP::Hep3Vector  _RightSOffsetInMu2e = _STMShieldingRef + CLHEP::Hep3Vector(R_dX, R_dY, R_dZ);

          stm._pSTMRightShieldingParams = std::unique_ptr<RightShielding>
          (new RightShielding(_RightShieldingBuild,
                              _RightS_Length,
                              _RightSleaddepth,
                              _RightScopperdepth,
                              _RightSBPdepth,
                              _RightSXmax,
                              _RightSOffsetInMu2e,
                              _RightSRotation));
    }


   ////////////////////////////////////////////////////////////////
   //STM Top Shielding

    if (_handstacked) {
      // The updated roof: three layers worked downward from the top of
      // the stack, already placed in parseConfig.
          stm._pTopWallParams = std::unique_ptr<TopWall>
          (new TopWall(_TopWallBuild,
                       _TopWallLeadLayers,
                       _TopWallSheets));
    } else {

      const double _TopS_Thickness = _TopSleaddepth*2 + _TopSBPdepth*2 + _TopScopperdepth;

      const double _T_Xlength1 = _TopSFront_LT + _TopSBarLeft + _TopSBarLeft + _TopSGapRight + _TopSGapRight;

      const double T_dX = -(_T_Xlength1 + _TopSXlength)/4 + _TopSBarLeft + _TopSGapLeft + _TopSFront_LT + (_STM_SSCdelta_WlL - _STM_SSCdelta_WlR)/2;
      const double T_dY = _FrontSHeightofRoom/2 +  _TopS_Thickness/2;
      const double T_dZ = _FrontS_Thickness + _TopSZlength/2 - _TopSZHole;

      const CLHEP::HepRotation _TopSRotation   = CLHEP::HepRotation::IDENTITY;
      const CLHEP::Hep3Vector  _TopSOffsetInMu2e = _STMShieldingRef + CLHEP::Hep3Vector(T_dX, T_dY, T_dZ);


          stm._pSTMTopShieldingParams = std::unique_ptr<TopShielding>
          (new TopShielding(_TopShieldingBuild,
                            _TopShieldingSkirtBuild,
                            _TopLiftBeam_L,
                            _TopLiftBeam_H,
                            _TopLiftBeam_T,
                            _TopLiftBeam_Xmove,
                            _TopSZlength,
                            _TopSXlength,
                            _TopSFront_LT,
                            _TopTFZlength,
                            _TopTFXlength,
                            _TopTBZlength,
                            _TopScontainerdepth,
                            _TopSleaddepth,
                            _TopScopperdepth,
                            _TopSBPdepth,
                            _TopSZHole,
                            _TopSBarLeft,
                            _TopSBarRight,
                            _TopSGapLeft,
                            _TopSGapRight,
                            _TopSLeak,
                            _TopSOffsetInMu2e,
                            _TopSRotation));
    }


   ////////////////////////////////////////////////////////////////
   //STM Inner Shielding

    if (_handstacked) {
      // The updated lining: prisms, boxes and bricks, each already
      // placed in parseConfig, and the bores they name.
          stm._pSTMInnerShieldingParams = std::unique_ptr<InnerShielding>
          (new InnerShielding(_InnerShieldingBuild,
                              _InnerShieldingPrisms,
                              _InnerShieldingBoxes,
                              _InnerShieldingBricks,
                              _InnerShieldingBores));
    } else {
      // The earlier description, whose geometry is written out in
      // constructSTM.cc, so only the switch is carried.
          stm._pSTMInnerShieldingParams = std::unique_ptr<InnerShielding>
          (new InnerShielding(_InnerShieldingBuild));
    }

   ////////////////////////////////////////////////////////////////
   //STM Back Shielding


      // Feeds only BackShielding's origin, which nothing reads (see
      // _STMShieldingRef above). _BottomSfloor_Zlength is set only for the
      // earlier description, so the updated one uses 0 rather than read it
      // unset.
      const double Back_dZ = _handstacked ? 0.
        : _BottomSfloor_Zlength + _BackSBPThick/2;

      const CLHEP::HepRotation _BackSRotation   = CLHEP::HepRotation::IDENTITY;
      const CLHEP::Hep3Vector  _BackSOffsetInMu2e = _STMShieldingRef + CLHEP::Hep3Vector(_BackS_dX, _BackS_dY, Back_dZ);

          stm._pSTMBackShieldingParams = std::unique_ptr<BackShielding>
          (new BackShielding(_BackShieldingBuild,
                              _BackSBPThick,
                              _BackSBPLength,
                              _BackSBPHeight,
                              _BackS_dX,
                              _BackS_dY,
                              _BackSPipeGap,
                              _BackSOffsetInMu2e,
                              _BackSRotation,
                              _BackSBPMaterial));

   ////////////////////////////////////////////////////////////////
   //STM Electronic Shielding

    if (_handstacked) {
      // Placed in the hall, not from the shield house: centered in x
      // xOffset from the beam axis, standing fromGround above the floor,
      // and set back from the east wall's inner face in z.
      std::vector<CLHEP::Hep3Vector> concreteCenters;
      std::vector<CLHEP::Hep3Vector> siGridCenters;
      if (_ElectronicShieldingBuild) {
        const double xCenter = BeamAxisAtEastWallInMu2e.x() + _ElectronicSXOffset;
        const double floorY  = FloorAtEastWallInMu2e.y();
        const double wallZ   = BeamAxisAtEastWallInMu2e.z();

        // The blocks stack upward from the floor, blockGap apart, their
        // downstream faces downstreamFaceToWall from the east wall.
        const double cdy = _ElectronicSConcreteDims.y();
        const double cdz = _ElectronicSConcreteDims.z();
        const double concreteBackZ = wallZ - _ElectronicSConcreteToWall;
        for (int i = 0; i < _ElectronicSConcreteBlockN; ++i) {
          concreteCenters.push_back(CLHEP::Hep3Vector(
            xCenter,
            floorY + _ElectronicSConcreteFromGround + i*(cdy + _ElectronicSConcreteBlockGap) + cdy/2,
            concreteBackZ - cdz/2));
        }

        // Each grid sits its gap behind the concrete's downstream face,
        // face to face, sharing the concrete's x center. The grids must
        // stay clear of each other and of the east wall.
        const double siY = floorY + _ElectronicSSiFromGround + _ElectronicSSiNY*_ElectronicSSiGridY/2;
        double prevBackZ = concreteBackZ;
        for (size_t i = 0; i < _ElectronicSSiGapsToConcrete.size(); ++i) {
          const double frontZ = concreteBackZ + _ElectronicSSiGapsToConcrete[i];
          if (frontZ < prevBackZ || frontZ + _ElectronicSSiGridZ > wallZ) {
            throw cet::exception("GEOM")
              << "STMMaker: stm.ElectronicShielding.si.gapsToConcrete entry "
              << i+1 << " (" << _ElectronicSSiGapsToConcrete[i] << " mm) puts"
              << " that silicon grid into the concrete, the grid before it, or"
              << " the east wall; the entries must increase by at least the"
              << " tile thickness (" << _ElectronicSSiGridZ << " mm) and stay"
              << " within " << _ElectronicSConcreteToWall - _ElectronicSSiGridZ
              << " mm.\n";
          }
          siGridCenters.push_back(CLHEP::Hep3Vector(xCenter, siY, frontZ + _ElectronicSSiGridZ/2));
          prevBackZ = frontZ + _ElectronicSSiGridZ;
        }

        // constructSTM places these in the downstream envelope, so each
        // must lie inside it; touching its floor is allowed.
        if (_stmDnStrEnvBuild) {
          const CLHEP::Hep3Vector envHalf(_stmDnStrEnvHalfWidth,
                                          _stmDnStrEnvHalfHeight,
                                          _stmDnStrEnvHalfLength);
          const CLHEP::Hep3Vector siHalf(_ElectronicSSiNX*_ElectronicSSiGridX/2,
                                         _ElectronicSSiNY*_ElectronicSSiGridY/2,
                                         _ElectronicSSiGridZ/2);
          auto checkInEnvelope = [&](CLHEP::Hep3Vector const & center,
                                     CLHEP::Hep3Vector const & half,
                                     std::string const & what) {
            // |d| + half <= envHalf on each axis is both faces inside.
            const CLHEP::Hep3Vector d = center - _stmDnStrEnvPositionInMu2e;
            const double tol = 1.e-6;
            if (std::abs(d.x()) + half.x() > envHalf.x() + tol ||
                std::abs(d.y()) + half.y() > envHalf.y() + tol ||
                std::abs(d.z()) + half.z() > envHalf.z() + tol) {
              throw cet::exception("GEOM")
                << "STMMaker: the electronic shielding's " << what
                << " (center " << center << ", half lengths " << half
                << ") does not fit in the STM downstream envelope (center "
                << _stmDnStrEnvPositionInMu2e << ", half lengths " << envHalf
                << ").\n";
            }
          };
          for (size_t i = 0; i < concreteCenters.size(); ++i) {
            checkInEnvelope(concreteCenters[i], _ElectronicSConcreteDims/2,
                            "concrete block " + std::to_string(i+1));
          }
          for (size_t i = 0; i < siGridCenters.size(); ++i) {
            checkInEnvelope(siGridCenters[i], siHalf,
                            "silicon grid " + std::to_string(i+1));
          }
        }

        if (_verbosityLevel > 0) {
          mf::LogInfo log("GEOM");
          log << "STMMaker: electronic shielding in Mu2e coordinates:\n";
          for (size_t i = 0; i < concreteCenters.size(); ++i) {
            log << "  concrete block " << i+1 << " center " << concreteCenters[i] << "\n";
          }
          for (size_t i = 0; i < siGridCenters.size(); ++i) {
            log << "  silicon grid " << i+1 << " center " << siGridCenters[i]
                << ", " << wallZ - siGridCenters[i].z() - _ElectronicSSiGridZ/2
                << " mm from the east wall\n";
          }
        }
      }

          stm._pSTMElectronicShieldingParams = std::unique_ptr<ElectronicShielding>
          (new ElectronicShielding(_ElectronicShieldingBuild,
                                   _ElectronicSConcreteDims/2,
                                   concreteCenters,
                                   CLHEP::Hep3Vector(_ElectronicSSiGridX,
                                                     _ElectronicSSiGridY,
                                                     _ElectronicSSiGridZ)/2,
                                   _ElectronicSSiNX,
                                   _ElectronicSSiNY,
                                   siGridCenters));
    } else {
          stm._pSTMElectronicShieldingParams = std::unique_ptr<ElectronicShielding>
          (new ElectronicShielding(_ElectronicShieldingBuild,
                                   _ElectronicSSiGridX,
                                   _ElectronicSSiGridY,
                                   _ElectronicSSiGridZ,
                                   _ElectronicSSiXcenter,
                                   _ElectronicSSiYcenter,
                                   _ElectronicSSiZcenter,
                                   _ElectronicSConcreteT,
                                   _ElectronicSGapToSi));
    }

   ////////////////////////////////////////////////////////////////
   //STM Absorber Shielding, or the SSC front shield that replaces it

    if (_handstacked) {
      // The updated geometry puts a stacked lead brick wall, a shelf and
      // two poly blocks where the earlier description had a single
      // absorber block.
      stm._pSSCFrontShieldParams = std::unique_ptr<SSCFrontShield>
              (new SSCFrontShield(_SSCFrontShieldBuild,
                                  _SSCFrontShieldBricks,
                                  _SSCFrontShieldShelfMaterial,
                                  _SSCFrontShieldShelfDim,
                                  _SSCFrontShieldShelfCenter,
                                  _SSCFrontShieldPoly1Material,
                                  _SSCFrontShieldPoly1Dim,
                                  _SSCFrontShieldPoly1Center,
                                  _SSCFrontShieldPoly1BoreR,
                                  _SSCFrontShieldPoly1BoreDX,
                                  _SSCFrontShieldPoly1BoreDY,
                                  _SSCFrontShieldPoly2Material,
                                  _SSCFrontShieldPoly2Dim,
                                  _SSCFrontShieldPoly2Center));
    } else {
          stm._pSTMSTM_AbsorberParams = std::unique_ptr<STM_Absorber>
          (new STM_Absorber(_STM_AbsorberBuild,
                            _STM_Absorber_hW,
                            _STM_Absorber_hH,
                            _STM_Absorber_hT,
                            _STM_Absorber_GaptoSSC));
    }

  }

  void STMMaker::parseConfig( SimpleConfig const & _config ){

    const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( _config, "stmMagnetField", "stm.magnet.field");


    _verbosityLevel            = _config.getInt("stm.verbosityLevel",0);

    // The updated, hand-stacked downstream shielding, taken from the NX STEP
    // export of the shield house. One switch for the whole downstream
    // geometry: those pieces are dimensioned against a shared datum and a
    // shared stack chain, so they ship together rather than one at a time.
    _handstacked               = _config.getBool("stm.downstream.handstacked",false);

    // The standard lead bricks. Read here rather than with any one
    // structure: most of the shield house is stacked from them, so they
    // are a shared primitive that many components refer to.
    if (_handstacked) {
      // Read by number, so adding a size to the geometry file is
      // enough: nothing here or downstream names a particular one.
      const int nBrickTypes = _config.getInt("stm.leadBrick.typeN");
      for (int t = 1; t <= nBrickTypes; ++t) {
        const std::string name =
          _config.getString("stm.leadBrick.type" + std::to_string(t));
        _leadBrickNames.push_back(name);
        _leadBrickDims.push_back(
          CLHEP::Hep3Vector(_config.getDouble("stm.leadBrick." + name + ".dx"),
                            _config.getDouble("stm.leadBrick." + name + ".dy"),
                            _config.getDouble("stm.leadBrick." + name + ".dz")));
      }
      _leadBrickWearY     = _config.getDouble("stm.leadBrick.wearY");
      _leadBrickWearXZ    = _config.getDouble("stm.leadBrick.wearXZ");
      _leadBrickMaterial  = _config.getString("stm.leadBrick.material");
    }

    _stmReferenceZ             = _config.getDouble("stm.referenceZ");  //was previously calculated automatically based on the location of the CRV-D

    // The sweeper magnet. Even when it is not built it is the anchor the
    // FOV collimator, the shield pipe and the magnet stand are placed
    // from, so halfLength and halfHeight are always required. The rest
    // is required only where something uses it: UpStrSpace when the
    // magnet is not placed off the pipe, the hole when the magnet or the
    // transport pipe through it is built, and the remainder only when
    // the magnet itself is built.
    _magnetBuild               = _config.getBool(  "stm.magnet.build",false);
    _magnetUpStrSpace          = readDoubleIfBuilt(_config, "stm.magnet.UpStrSpace",
                                                   _magnetBuild || !_config.getBool("stm.magnet.usePipeAsOrigin", false));
    _magnetHalfLength          = _config.getDouble("stm.magnet.halfLength");
    _magnetHalfWidth           = readDoubleIfBuilt(_config, "stm.magnet.halfWidth", _magnetBuild);
    _magnetHalfHeight          = _config.getDouble("stm.magnet.halfHeight");
    const bool magnetHoleNeeded = _magnetBuild || _config.getBool("stm.pipe.build", false);
    _magnetHoleHalfWidth       = readDoubleIfBuilt(_config, "stm.magnet.holeHalfWidth",  magnetHoleNeeded);
    _magnetHoleHalfHeight      = readDoubleIfBuilt(_config, "stm.magnet.holeHalfHeight", magnetHoleNeeded);
    _magnetHoleXOffset         = _config.getDouble("stm.magnet.holeXOffset", 0.);
    _magnetHoleYOffset         = _config.getDouble("stm.magnet.holeYOffset", 0.);
    _magnetMaterial            = readStringIfBuilt(_config, "stm.magnet.material", _magnetBuild);
    _magnetHasLiner            = _config.getBool("stm.magnet.hasLiner", true);
    _magnetField               = readDoubleIfBuilt(_config, "stm.magnet.field", _magnetBuild);
    //_magnetFieldVisible        = _config.getBool(  "stm.magnet.fieldVisible",false);
    _magnetFieldVisible        = geomOptions->isVisible("stmMagnetField");

    _FOVCollimatorBuild            = _config.getBool(  "stm.FOVcollimator.build");
    _FOVCollimatorMaterial         = _config.getString("stm.FOVcollimator.material");
    _FOVCollimatorUpStrSpace       = _config.getDouble("stm.FOVcollimator.UpStrSpace");
    _FOVCollimatorHalfWidth        = _config.getDouble("stm.FOVcollimator.halfWidth");
    _FOVCollimatorHalfHeight       = _config.getDouble("stm.FOVcollimator.halfHeight");
    _FOVCollimatorHalfLength       = _config.getDouble("stm.FOVcollimator.halfLength");
    _FOVCollimatorLinerBuild       = _config.getBool(  "stm.FOVcollimator.liner.build");
    // The liner's size and the cutout it sits in also place the poly
    // absorber, so they are required when either is built. The absorber
    // switch is read here with the same fallback constructSTM gives it.
    const bool FOVLinerNeeded = _FOVCollimatorLinerBuild ||
                                _config.getBool("stm.FOVcollimator.absorber.build", false);
    _FOVCollimatorLinerMaterial    = readStringIfBuilt(_config, "stm.FOVcollimator.liner.material", _FOVCollimatorLinerBuild);
    _FOVCollimatorLinerHalfWidth   = readDoubleIfBuilt(_config, "stm.FOVcollimator.liner.halfWidth", FOVLinerNeeded);
    _FOVCollimatorLinerHalfHeight  = readDoubleIfBuilt(_config, "stm.FOVcollimator.liner.halfHeight", FOVLinerNeeded);
    _FOVCollimatorLinerHalfLength  = readDoubleIfBuilt(_config, "stm.FOVcollimator.liner.halfLength", FOVLinerNeeded);
    _FOVCollimatorLinerCutOutHalfLength  = readDoubleIfBuilt(_config, "stm.FOVcollimator.liner.cutOutHalfLength", FOVLinerNeeded);
    _FOVCollimatorHole1RadiusUpStr = _config.getDouble("stm.FOVcollimator.hole1.radiusUpStr");
    // The FOV collimator is a box with one centered hole of radius
    // hole1.radiusUpStr. STMCollimator is shared with the earlier SS
    // collimator and carries more, but nothing asks the FOV collimator
    // for the rest, so it is not read from the geometry file.
    _FOVCollimatorHole1xOffset        = 0.;
    _FOVCollimatorHole1RadiusDnStr    = 0.;
    _FOVCollimatorHole1LinerBuild     = false;
    _FOVCollimatorHole1LinerThickness = 0.;
    _FOVCollimatorHole2Build          = false;
    _FOVCollimatorHole2xOffset        = 0.;
    _FOVCollimatorHole2RadiusUpStr    = 0.;
    _FOVCollimatorHole2RadiusDnStr    = 0.;
    _FOVCollimatorHole2LinerBuild     = false;
    _FOVCollimatorHole2LinerThickness = 0.;
    _FOVCollimatorHoleLinerMaterial   = "";

    // The transport pipe through the magnet. Its keys are required only
    // when it is built; nothing else is placed from it.
    _pipeBuild                 = _config.getBool(  "stm.pipe.build",false);
    _pipeRadiusIn              = readDoubleIfBuilt(_config, "stm.pipe.rIn", _pipeBuild);
    _pipeRadiusOut             = readDoubleIfBuilt(_config, "stm.pipe.rOut", _pipeBuild);
    _pipeMaterial              = readStringIfBuilt(_config, "stm.pipe.material", _pipeBuild);
    _pipeGasMaterial           = readStringIfBuilt(_config, "stm.pipe.gas.material", _pipeBuild);
    _pipeUpStrSpace            = readDoubleIfBuilt(_config, "stm.pipe.UpStrSpace", _pipeBuild);
    _pipeDnStrHalfLength       = readDoubleIfBuilt(_config, "stm.pipe.DnStrHalfLength", _pipeBuild);
    _pipeUpStrWindowMaterial   = readStringIfBuilt(_config, "stm.pipe.UpStrWindow.material", _pipeBuild);
    _pipeUpStrWindowHalfLength = readDoubleIfBuilt(_config, "stm.pipe.UpStrWindow.halfLength", _pipeBuild);
    _pipeDnStrWindowMaterial   = readStringIfBuilt(_config, "stm.pipe.DnStrWindow.material", _pipeBuild);
    _pipeDnStrWindowHalfLength = readDoubleIfBuilt(_config, "stm.pipe.DnStrWindow.halfLength", _pipeBuild);
    _pipeFlangeHalfLength      = readDoubleIfBuilt(_config, "stm.pipe.flange.halfLength", _pipeBuild);
    _pipeFlangeOverhangR       = readDoubleIfBuilt(_config, "stm.pipe.flange.overhangR", _pipeBuild);

    _magnetTableBuild          = _config.getBool(  "stm.magnet.stand.build",false);
    _magnetTableMaterial       = _config.getString("stm.magnet.stand.material");
    _magnetTableTopExtraWidth  = _config.getDouble("stm.magnet.stand.topExtraWidth");
    _magnetTableTopExtraLength = _config.getDouble("stm.magnet.stand.topExtraLength");
    _magnetTableTopHalfHeight  = _config.getDouble("stm.magnet.stand.topHalfHeight");
    _magnetTableLegRadius      = _config.getDouble("stm.magnet.stand.legRadius");

    // The earlier lead SS collimator, its detector stand and the two Ge
    // detectors, none of which the current geometry builds. Their keys
    // are required only when the component is built.
    _SSCollimatorBuild            = _config.getBool(  "stm.SScollimator.build",false);
    _SSCollimatorMaterial         = readStringIfBuilt(_config, "stm.SScollimator.material", _SSCollimatorBuild);
    _SSCollimatorUpStrSpace       = readDoubleIfBuilt(_config, "stm.SScollimator.UpStrSpace", _SSCollimatorBuild);
    _SSCollimatorHalfWidth        = readDoubleIfBuilt(_config, "stm.SScollimator.halfWidth", _SSCollimatorBuild);
    _SSCollimatorHalfHeight       = readDoubleIfBuilt(_config, "stm.SScollimator.halfHeight", _SSCollimatorBuild);
    _SSCollimatorHalfLength       = readDoubleIfBuilt(_config, "stm.SScollimator.halfLength", _SSCollimatorBuild);
    _SSCollimatorLinerBuild       = readBoolIfBuilt(_config, "stm.SScollimator.liner.build", _SSCollimatorBuild);
    _SSCollimatorLinerMaterial    = readStringIfBuilt(_config, "stm.SScollimator.liner.material", _SSCollimatorBuild);
    _SSCollimatorLinerHalfWidth   = readDoubleIfBuilt(_config, "stm.SScollimator.liner.halfWidth", _SSCollimatorBuild);
    _SSCollimatorLinerHalfHeight  = readDoubleIfBuilt(_config, "stm.SScollimator.liner.halfHeight", _SSCollimatorBuild);
    _SSCollimatorLinerHalfLength  = readDoubleIfBuilt(_config, "stm.SScollimator.liner.halfLength", _SSCollimatorBuild);
    _SSCollimatorLinerCutOutHalfLength  = readDoubleIfBuilt(_config, "stm.SScollimator.liner.cutOutHalfLength", _SSCollimatorBuild);
    _SSCollimatorHole1xOffset     = readDoubleIfBuilt(_config, "stm.SScollimator.hole1.xoffset", _SSCollimatorBuild);
    _SSCollimatorHole1RadiusUpStr = readDoubleIfBuilt(_config, "stm.SScollimator.hole1.radiusUpStr", _SSCollimatorBuild);
    _SSCollimatorHole1RadiusDnStr = readDoubleIfBuilt(_config, "stm.SScollimator.hole1.radiusDnStr", _SSCollimatorBuild);
    _SSCollimatorHole1LinerBuild     = readBoolIfBuilt(_config, "stm.SScollimator.hole1.liner.build", _SSCollimatorBuild);
    _SSCollimatorHole1LinerThickness = readDoubleIfBuilt(_config, "stm.SScollimator.hole1.liner.thickness", _SSCollimatorBuild);
    _SSCollimatorHole2Build       = readBoolIfBuilt(_config, "stm.SScollimator.hole2.build", _SSCollimatorBuild);
    _SSCollimatorHole2xOffset     = readDoubleIfBuilt(_config, "stm.SScollimator.hole2.xoffset", _SSCollimatorBuild);
    _SSCollimatorHole2RadiusUpStr = readDoubleIfBuilt(_config, "stm.SScollimator.hole2.radiusUpStr", _SSCollimatorBuild);
    _SSCollimatorHole2RadiusDnStr = readDoubleIfBuilt(_config, "stm.SScollimator.hole2.radiusDnStr", _SSCollimatorBuild);
    _SSCollimatorHole2LinerBuild     = readBoolIfBuilt(_config, "stm.SScollimator.hole2.liner.build", _SSCollimatorBuild);
    _SSCollimatorHole2LinerThickness = readDoubleIfBuilt(_config, "stm.SScollimator.hole2.liner.thickness", _SSCollimatorBuild);
    _SSCollimatorHoleLinerMaterial = readStringIfBuilt(_config, "stm.SScollimator.hole.liner.material", _SSCollimatorBuild);

    _detectorTableBuild          = _config.getBool(  "stm.detector.stand.build",false);
    _detectorTableMaterial       = readStringIfBuilt(_config, "stm.detector.stand.material", _detectorTableBuild);
    _detectorTableTopExtraWidth  = readDoubleIfBuilt(_config, "stm.detector.stand.topExtraWidth", _detectorTableBuild);
    _detectorTableTopExtraLength = readDoubleIfBuilt(_config, "stm.detector.stand.topExtraLength", _detectorTableBuild);
    _detectorTableTopHalfHeight  = readDoubleIfBuilt(_config, "stm.detector.stand.topHalfHeight", _detectorTableBuild);
    _detectorTableLegRadius      = readDoubleIfBuilt(_config, "stm.detector.stand.legRadius", _detectorTableBuild);

    _detector1Build                    = _config.getBool(  "stm.det1.build",false);
    _detector1CrystalMaterial          = readStringIfBuilt(_config, "stm.det1.material", _detector1Build);
    _detector1CrystalRadiusIn          = readDoubleIfBuilt(_config, "stm.det1.rIn", _detector1Build);
    _detector1CrystalRadiusOut         = readDoubleIfBuilt(_config, "stm.det1.rOut", _detector1Build);
    _detector1CrystalHalfLength        = readDoubleIfBuilt(_config, "stm.det1.halfLength", _detector1Build);
    _detector1xOffset                  = readDoubleIfBuilt(_config, "stm.det1.xoffset", _detector1Build);
    _detector1CanMaterial              = readStringIfBuilt(_config, "stm.det1.can.material", _detector1Build);
    _detector1CanRadiusIn              = readDoubleIfBuilt(_config, "stm.det1.can.rIn", _detector1Build);
    _detector1CanRadiusOut             = readDoubleIfBuilt(_config, "stm.det1.can.rOut", _detector1Build);
    _detector1CanHalfLength            = readDoubleIfBuilt(_config, "stm.det1.can.halfLength", _detector1Build);
    _detector1CanUpStrSpace            = readDoubleIfBuilt(_config, "stm.det1.can.UpStrSpace", _detector1Build);
    _detector1CanUpStrWindowMaterial   = readStringIfBuilt(_config, "stm.det1.can.UpStrWindowMaterial", _detector1Build);
    _detector1CanUpStrWindowHalfLength = readDoubleIfBuilt(_config, "stm.det1.can.UpStrWindowHalfLength", _detector1Build);
    _detector1CanGasMaterial           = readStringIfBuilt(_config, "stm.det1.can.gas", _detector1Build);

    _detector2Build                    = _config.getBool(  "stm.det2.build",false);
    _detector2CrystalMaterial          = readStringIfBuilt(_config, "stm.det2.material", _detector2Build);
    _detector2CrystalRadiusIn          = readDoubleIfBuilt(_config, "stm.det2.rIn", _detector2Build);
    _detector2CrystalRadiusOut         = readDoubleIfBuilt(_config, "stm.det2.rOut", _detector2Build);
    _detector2CrystalHalfLength        = readDoubleIfBuilt(_config, "stm.det2.halfLength", _detector2Build);
    _detector2xOffset                  = readDoubleIfBuilt(_config, "stm.det2.xoffset", _detector2Build);
    _detector2CanMaterial              = readStringIfBuilt(_config, "stm.det2.can.material", _detector2Build);
    _detector2CanRadiusIn              = readDoubleIfBuilt(_config, "stm.det2.can.rIn", _detector2Build);
    _detector2CanRadiusOut             = readDoubleIfBuilt(_config, "stm.det2.can.rOut", _detector2Build);
    _detector2CanHalfLength            = readDoubleIfBuilt(_config, "stm.det2.can.halfLength", _detector2Build);
    _detector2CanUpStrSpace            = readDoubleIfBuilt(_config, "stm.det2.can.UpStrSpace", _detector2Build);
    _detector2CanUpStrWindowMaterial   = readStringIfBuilt(_config, "stm.det2.can.UpStrWindowMaterial", _detector2Build);
    _detector2CanUpStrWindowHalfLength = readDoubleIfBuilt(_config, "stm.det2.can.UpStrWindowHalfLength", _detector2Build);
    _detector2CanGasMaterial           = readStringIfBuilt(_config, "stm.det2.can.gas", _detector2Build);

    // The z budget measured upstream from the east hall wall. Only the
    // four components above are placed from it.
    _stmZAllowed = readDoubleIfBuilt(_config, "stm.z.allowed",
                                     _SSCollimatorBuild || _detectorTableBuild ||
                                     _detector1Build    || _detector2Build);

    _shieldBuild                = _config.getBool(  "stm.shield.build",false);
    _shieldRadiusIn             = _config.getDouble("stm.shield.rIn");
    _shieldHasLiner             = _config.getBool(  "stm.shield.hasLiner", true /*true default for backwards compatibility*/);
    _shieldLinerWidth           = _config.getDouble("stm.shield.widthLiner");
    _shieldRadiusOut            = _config.getDouble("stm.shield.rOut");
    _shieldPipeHalfLength       = _config.getDouble("stm.shield.pipe.halfLength");
    _shieldMaterialLiner        = _config.getString("stm.shield.materialLiner");
    _shieldMaterial             = _config.getString("stm.shield.material");
    _shieldMatchPipeBlock       = _config.getBool  ("stm.shield.matchPipeBlock", false);
    _shieldUpStrSpace           = _config.getDouble("stm.shield.UpStrSpace");
    _shieldDnStrSpace           = _config.getDouble("stm.shield.DnStrSpace");
    // The mating block (DnStrWall). Its thickness is required only when
    // it is built; without it the block has no extent, so it adds no
    // length to the magnet stand and no offset to the magnet when the
    // pipe is not matched to it. The other DnStrWall keys already fall
    // back when absent. The switch itself still defaults to true for
    // older geometry files, so a file without the block must say so.
    _shieldBuildMatingBlock     = _config.getBool("stm.shield.matingBlock.build", true); // default to true because that is what older versions did
    _shieldDnStrWallHalfLength  = readDoubleIfBuilt(_config, "stm.shield.DnStrWall.halfLength", _shieldBuildMatingBlock);
    _shieldDnStrWallHoleRadius  = _config.getDouble("stm.shield.DnStrWall.holeRadius", -1.);
    _shieldDnStrWallHalfHeight  = _config.getDouble("stm.shield.DnStrWall.halfHeight", -1.);
    _shieldDnStrWallHalfWidth   = _config.getDouble("stm.shield.DnStrWall.halfWidth", -1.);
    _shieldDnStrWallGap         = _config.getDouble("stm.shield.DnStrWall.gap", 0.);
    _shieldUpStrWallGap         = _config.getDouble("stm.shield.UpStrWall.gap", 0.); //only if using pipe as origin
    _shieldDnStrWallMaterial    = _config.getString("stm.shield.DnStrWall.material", _shieldMaterial);
    _shieldPipeUpStrAirGap      = _config.getDouble("stm.shield.pipe.upStrAirGap", 0); // default to 0 for backwards compatibility

    _stmDnStrEnvBuild       = _config.getBool("stm.downstream.build");
    _stmDnStrEnvHalfLength  = _config.getDouble("stm.downstream.halfLength");
    _stmDnStrEnvHalfWidth   = _config.getDouble("stm.downstream.halfWidth");
    _stmDnStrEnvHalfHeight  = - _config.getDouble("yOfFloorSurface.below.mu2eOrigin");
    _stmDnStrEnvMaterial    = _config.getString("stm.downstream.material");

    _STM_SSCBuild          = _config.getBool(  "stm.STM_SSC.build");
    _STM_SSCVDBuild        = _config.getBool(  "stm.STM_SSC.VDbuild");
    _STM_SSCW_height       = _config.getDouble("stm.STM_SSC.W_height");
    _STM_SSCWdepth_f       = _config.getDouble("stm.STM_SSC.Wdepth_f");
    _STM_SSCWdepth_b       = _config.getDouble("stm.STM_SSC.Wdepth_b");

    // Each description reads only the keys it has: the earlier SSC is a
    // middle block plus two wings with the bores given as aperture areas,
    // the updated one a single block with the bores given as radii. Reading
    // them apart means a geometry file need not carry the other's parameters.
    if (_handstacked) {
      _STM_SSCW_width      = _config.getDouble("stm.STM_SSC.W_width");
      _STM_SSCr_LaBr_f     = _config.getDouble("stm.STM_SSC.r_LaBr_f");
      _STM_SSCr_HPGe_f     = _config.getDouble("stm.STM_SSC.r_HPGe_f");
      _STM_SSCr_LaBr_b     = _config.getDouble("stm.STM_SSC.r_LaBr_b");
      _STM_SSCr_HPGe_b     = _config.getDouble("stm.STM_SSC.r_HPGe_b");
      // The height of the bore above the baseplate: the datum every
      // hand-stacked piece downstream measures its own height from.
      _STM_SSCboreToBase   = _config.getDouble("stm.STM_SSC.boreToBase");
    } else {
      _STM_SSCdelta_WlR      = _config.getDouble("stm.STM_SSC.delta_WlR");
      _STM_SSCdelta_WlL      = _config.getDouble("stm.STM_SSC.delta_WlL");
      _STM_SSCW_middle       = _config.getDouble("stm.STM_SSC.W_middle");
      _STM_SSCAperture_HPGe1 = _config.getDouble("stm.STM_SSC.Aperture_HPGe1");
      _STM_SSCAperture_HPGe2 = _config.getDouble("stm.STM_SSC.Aperture_HPGe2");
      _STM_SSCAperture_LaBr1 = _config.getDouble("stm.STM_SSC.Aperture_LaBr1");
      _STM_SSCAperture_LaBr2 = _config.getDouble("stm.STM_SSC.Aperture_LaBr2");
    }
    _STM_SSCoffset_Spot    = _config.getDouble("stm.STM_SSC.offset_Spot");
    _STM_SSCleak           = _config.getDouble("stm.STM_SSC.leak");
    _STM_SSCFrontToWall    = _config.getDouble("stm.STM_SSC.FrontToWall");
    _STM_SSCZGap           = _config.getDouble("stm.STM_SSC.ZGap");
    _STM_SSCZGapBack       = _config.getDouble("stm.STM_SSC.ZGapBack");
    _STM_SSCMaterial       = _config.getString("stm.STM_SSC.material");

    _SSCSupportBuild          = _config.getBool(  "stm.SSCSupport.build");

    // The updated support is a five-plate steel cradle, not the earlier
    // table/legs/base/walls plus shims, so the two read different keys.
    if (_handstacked) {
      _SSCSupportdepth        = _config.getDouble("stm.SSCSupport.depth");
      _SSCSupportside_T       = _config.getDouble("stm.SSCSupport.side_T");
      _SSCSupportside_H       = _config.getDouble("stm.SSCSupport.side_H");
      _SSCSupportplate_base_T = _config.getDouble("stm.SSCSupport.base_T");
      _SSCSupportbottom_T     = _config.getDouble("stm.SSCSupport.bottom_T");
      _SSCSupporttop_T        = _config.getDouble("stm.SSCSupport.top_T");
      _SSCSupportMaterial     = _config.getString("stm.SSCSupport.material");
    } else {
      _SSCSupporttable_L      = _config.getDouble("stm.SSCSupport.table_L");
      _SSCSupporttable_H      = _config.getDouble("stm.SSCSupport.table_H");
      _SSCSupporttable_T      = _config.getDouble("stm.SSCSupport.table_T");
      _SSCSupportleg_L      = _config.getDouble("stm.SSCSupport.leg_L");
      _SSCSupportleg_H      = _config.getDouble("stm.SSCSupport.leg_H");
      _SSCSupportleg_T      = _config.getDouble("stm.SSCSupport.leg_T");
      _SSCSupportbase_L      = _config.getDouble("stm.SSCSupport.base_L");
      _SSCSupportbase_H      = _config.getDouble("stm.SSCSupport.base_H");
      _SSCSupportbase_T      = _config.getDouble("stm.SSCSupport.base_T");
      _SSCSupportwall_L      = _config.getDouble("stm.SSCSupport.wall_L");
      _SSCSupportwall_H      = _config.getDouble("stm.SSCSupport.wall_H");
      _SSCSupportwall_T      = _config.getDouble("stm.SSCSupport.wall_T");
      _SSCSupporthole_H      = _config.getDouble("stm.SSCSupport.hole_H");
      _SSCSupporthole_T      = _config.getDouble("stm.SSCSupport.hole_T");
      _SSCSupportFLeadStand_L      = _config.getDouble("stm.SSCSupport.FLeadStand_L");
      _SSCSupportFLeadStand_H      = _config.getDouble("stm.SSCSupport.FLeadStand_H");
      _SSCSupportFLeadStand_T      = _config.getDouble("stm.SSCSupport.FLeadStand_T");
      _SSCSupportFLeadShim_H       = _config.getDouble("stm.SSCSupport.FLeadShim_H");
      _SSCSupportFLeadShim_T       = _config.getDouble("stm.SSCSupport.FLeadShim_T");
      _SSCSupportFAluminumShim_T   = _config.getDouble("stm.SSCSupport.FAluminumShim_T");
      _SSCSupportFAluminumExtra_L  = _config.getDouble("stm.SSCSupport.FAluminumExtra_L");
      _SSCSupportFAluminumExtra_H  = _config.getDouble("stm.SSCSupport.FAluminumExtra_H");
    }

    // ---- SSC front shield ---------------------------------------------
    //
    // The geometry file writes every center as if the structure sat
    // exactly on the SSC axis. Three offsets displace it, and they are
    // resolved here so that nothing downstream has to know which one
    // applies to which piece:
    //
    //   offsetX/offsetY        the whole structure
    //   LaBrSideBrickOffsetX   the bricks, except those a group's
    //                          followPoly2 names, which are constrained
    //                          by poly2 instead
    //   poly2.offsetX          poly2, and those same bricks
    if (_handstacked) {

      _SSCFrontShieldBuild = _config.getBool("stm.SSCFrontShield.build");

      const double fsOffX      = _config.getDouble("stm.SSCFrontShield.offsetX");
      const double fsOffY      = _config.getDouble("stm.SSCFrontShield.offsetY");
      const double fsLaBrOffX  = _config.getDouble("stm.SSCFrontShield.LaBrSideBrickOffsetX");
      const double fsPoly2OffX = _config.getDouble("stm.SSCFrontShield.poly2.offsetX");

      const CLHEP::Hep3Vector fsShift(fsOffX, fsOffY, 0.);

      // Read one group of bricks: a type, n centers, and either one
      // orientation for the whole group or one per brick. extraX is the
      // x offset each brick of the group follows, on top of the
      // whole-structure shift.
      auto readBrickGroup = [&](std::string const & base, int n,
                                std::vector<double> const & extraX) {
        const int type = _config.getInt(base + ".type");
        if (type < 1 || type > int(_leadBrickDims.size())) {
          throw cet::exception("GEOM")
            << "STMMaker: " << base << ".type names brick type " << type
            << ", but stm.leadBrick.typeN defines only "
            << _leadBrickDims.size() << ".\n";
        }
        const std::string shared = _config.getString(base + ".orientation", "");

        std::vector<BrickWallBrick> out;
        out.reserve(n);
        for (int i = 1; i <= n; ++i) {
          std::ostringstream ckey, okey;
          ckey << base << ".center"      << i;
          okey << base << ".orientation" << i;

          std::vector<double> c;
          _config.getVectorDouble(ckey.str(), c, 3);

          BrickWallBrick b;
          b.type        = type;
          b.orientation = shared.empty()
                        ? _config.getString(okey.str(), "000") : shared;
          b.center      = CLHEP::Hep3Vector(c[0] + extraX.at(i-1), c[1], c[2])
                        + fsShift;
          out.push_back(b);
        }
        return out;
      };

      _SSCFrontShieldBricks.clear();

      const int fsGroupN = _config.getInt("stm.SSCFrontShield.brickGroupN");
      for (int g = 1; g <= fsGroupN; ++g) {
        std::ostringstream base;
        base << "stm.SSCFrontShield.brickGroup" << g;
        const int n = _config.getInt(base.str() + ".n");

        // Most bricks follow the LaBr side stack; the few that sit
        // against poly2 follow its offset instead, and the geometry
        // file names them by position rather than this code assuming
        // which group they are in.
        std::vector<int> followPoly2;
        _config.getVectorInt(base.str() + ".followPoly2", followPoly2,
                             std::vector<int>());

        std::vector<double> extraX(n, fsLaBrOffX);
        for (int i : followPoly2) {
          if (i < 1 || i > n) {
            throw cet::exception("GEOM")
              << "STMMaker: " << base.str() << ".followPoly2 names brick "
              << i << ", but the group has " << n << ".\n";
          }
          extraX[i-1] = fsPoly2OffX;
        }

        const std::vector<BrickWallBrick> group =
          readBrickGroup(base.str(), n, extraX);
        _SSCFrontShieldBricks.insert(_SSCFrontShieldBricks.end(),
                                     group.begin(), group.end());
      }

      // The shelf and poly1 take only the whole-structure offset.
      _SSCFrontShieldShelfMaterial = _config.getString("stm.SSCFrontShield.shelf.material");
      _SSCFrontShieldShelfDim = CLHEP::Hep3Vector(_config.getDouble("stm.SSCFrontShield.shelf.dx"),
                                                  _config.getDouble("stm.SSCFrontShield.shelf.dy"),
                                                  _config.getDouble("stm.SSCFrontShield.shelf.dz"));
      {
        std::vector<double> c;
        _config.getVectorDouble("stm.SSCFrontShield.shelf.center", c, 3);
        _SSCFrontShieldShelfCenter = CLHEP::Hep3Vector(c[0], c[1], c[2]) + fsShift;
      }

      _SSCFrontShieldPoly1Material = _config.getString("stm.SSCFrontShield.poly1.material");
      _SSCFrontShieldPoly1Dim = CLHEP::Hep3Vector(_config.getDouble("stm.SSCFrontShield.poly1.dx"),
                                                  _config.getDouble("stm.SSCFrontShield.poly1.dy"),
                                                  _config.getDouble("stm.SSCFrontShield.poly1.dz"));
      {
        std::vector<double> c;
        _config.getVectorDouble("stm.SSCFrontShield.poly1.center", c, 3);
        _SSCFrontShieldPoly1Center = CLHEP::Hep3Vector(c[0], c[1], c[2]) + fsShift;
      }
      _SSCFrontShieldPoly1BoreR  = _config.getDouble("stm.SSCFrontShield.poly1.boreR");
      _SSCFrontShieldPoly1BoreDX = _config.getDouble("stm.SSCFrontShield.poly1.boreDX");
      _SSCFrontShieldPoly1BoreDY = _config.getDouble("stm.SSCFrontShield.poly1.boreDY");

      // poly2 carries its own x offset on top of the structure offset.
      _SSCFrontShieldPoly2Material = _config.getString("stm.SSCFrontShield.poly2.material");
      _SSCFrontShieldPoly2Dim = CLHEP::Hep3Vector(_config.getDouble("stm.SSCFrontShield.poly2.dx"),
                                                  _config.getDouble("stm.SSCFrontShield.poly2.dy"),
                                                  _config.getDouble("stm.SSCFrontShield.poly2.dz"));
      {
        std::vector<double> c;
        _config.getVectorDouble("stm.SSCFrontShield.poly2.center", c, 3);
        _SSCFrontShieldPoly2Center = CLHEP::Hep3Vector(c[0] + fsPoly2OffX, c[1], c[2]) + fsShift;
      }
    }


    if (_handstacked) {
      // ---- Front shielding, right ---------------------------------------
      //
      // The config gives a pattern, not positions: courses as lists of
      // brick types, a layer sequence in depth, and a map saying which
      // bores pass through which brick. Everything is expanded here, so
      // the construction code receives finished placements.
      //
      // Three constraints fix the section. brickEndX pins the +x end of
      // every course; the bore sits boreToBase above the baseplate with
      // courses on a 4 in pitch from it; and wallToCradleGap sets the
      // wall's front face behind the cradle. Pieces then butt together.

      _FrontShieldingRightBuild = _config.getBool("stm.FrontShieldingRight.build");

      const double fsrBrickEndX = _config.getDouble("stm.FrontShieldingRight.brickEndX");
      const double fsrGap       = _config.getDouble("stm.FrontShieldingRight.wallToCradleGap");
      const double fsrPitch     = _config.getDouble("stm.FrontShieldingRight.coursePitch");
      const std::string fsrBrickOrient =
        _config.getString("stm.FrontShieldingRight.brickOrientation");

      // The wall starts behind the cradle and grows downstream.
      const double fsrFrontZ = _SSCSupportdepth + fsrGap;

      // Kept for the sections that butt this one; see STMMaker.hh.
      _FrontShieldingRightBrickEndX = fsrBrickEndX;
      _FrontShieldingRightFrontZ    = fsrFrontZ;

      const std::vector<BrickWallBore> fsrBores =
        readBores(_config, "stm.FrontShieldingRight", _STM_SSCoffset_Spot);
      const int fsrBoreN = int(fsrBores.size());

      // Which bores go through which brick: {layer, course, position,
      // bore}, one entry per hole. `used` is set when the entry matches a
      // brick as the layers are laid; an entry naming a layer, course or
      // position that has no brick would otherwise be skipped silently,
      // leaving that brick solid.
      struct FSRBoreMapEntry { int layer, course, position, bore; bool used; };
      std::vector<FSRBoreMapEntry> fsrBoreMap;
      const int fsrBoreMapN = _config.getInt("stm.FrontShieldingRight.leadBoreMapN");
      for (int i = 1; i <= fsrBoreMapN; ++i) {
        std::ostringstream key;
        key << "stm.FrontShieldingRight.leadBoreMap" << i;
        std::vector<int> e;
        _config.getVectorInt(key.str(), e, 4);
        fsrBoreMap.push_back({e[0], e[1], e[2], e[3], false});
      }

      // Courses, numbered so that n = 0 is the one the beam passes
      // through.
      std::vector<int> fsrCourse;
      const int fsrCourseN = _config.getInt("stm.FrontShieldingRight.courseN");
      for (int i = 1; i <= fsrCourseN; ++i) {
        std::ostringstream key;
        key << "stm.FrontShieldingRight.course" << i;
        fsrCourse.push_back(_config.getInt(key.str()));
      }

      // Where the courses actually sit.
      //
      // They stack up from the baseplate -- block, then course after
      // course -- so their height comes from that chain, not from the
      // beam. The beam's own height is boreToBase above the same
      // baseplate, and the two coincide only if the chain happens to
      // leave half a course on the beam plane. That is true of this
      // section as built (236.000 - 83.600 - 101.600 = 50.800, half a
      // course) but it is arithmetic, not a constraint, so the offset
      // is derived here and checked against the bores below.
      const int fsrCoursesBelow =
        std::count_if(fsrCourse.begin(), fsrCourse.end(),
                      [](int n) { return n < 0; });
      const double fsrBlockDyForCourses =
        _config.getDouble("stm.FrontShieldingRight.block.dy");
      const double fsrCourseOffsetY =
        -_STM_SSCboreToBase + fsrBlockDyForCourses
        + fsrCoursesBelow*fsrPitch + fsrPitch/2;

      // How much of the course a brick of this type takes up. Rotation
      // is also taken into consideration
      const CLHEP::Hep3Vector fsrCourseDir(-1., 0., 0.);  // courses run toward -x
      const CLHEP::Hep3Vector fsrPitchDir  ( 0., 1., 0.);  // courses stack upward
      const CLHEP::Hep3Vector fsrDepthDir  ( 0., 0., 1.);  // layers face downstream
      CLHEP::HepRotation fsrBrickRot(CLHEP::HepRotation::IDENTITY);
      {
        OrientationResolver OR;
        OR.getRotationFromOrientation(fsrBrickRot, fsrBrickOrient);
      }
      auto fsrTypeWidth = [&](int type) {
        return brickSpan(_leadBrickDims, type, fsrBrickRot, fsrCourseDir,
                         "stm.FrontShieldingRight");
      };

      // How much wear comes off each face along a given direction. The
      // wear is not isotropic -- one value for whichever of the
      // brick's own axes ends up vertical, another for the other two --
      // so it has to be projected through the rotation just as the
      // dimensions are. Needed by the bore guard, which has to compare
      // against the solid the hole is actually cut from, not the
      // as-delivered brick. The brick's own axis k lands on
      // rot^-1 * e_k in Mu2e, as in spanAlong.
      auto fsrWearAlong = [&](CLHEP::Hep3Vector const & dir) {
        double w = 0.;
        for (int k = 0; k < 3; ++k) {
          const CLHEP::Hep3Vector axis =
            fsrBrickRot.inverse() * CLHEP::Hep3Vector(k == 0, k == 1, k == 2);
          w += std::abs(axis.dot(dir))
             * ((std::abs(axis.y()) > 0.5) ? _leadBrickWearY : _leadBrickWearXZ);
        }
        return w;
      };

      // Walk the layer sequence in depth. A lead layer expands into a
      // BrickWall; a sheet layer into plates held by the section.
      const int fsrLayerN = _config.getInt("stm.FrontShieldingRight.layerN");
      double fsrZ = fsrFrontZ;
      int fsrLeadLayer = 0;

      const std::string fsrSheetMat =
        _config.getString("stm.FrontShieldingRight.BPlayer.material");
      const double fsrSheetT  = _config.getDouble("stm.FrontShieldingRight.BPlayer.thickness");
      const double fsrSheetLo = _config.getDouble("stm.FrontShieldingRight.BPlayer.dyLower");
      const double fsrSheetUp = _config.getDouble("stm.FrontShieldingRight.BPlayer.dyUpper");

      // The courses span this, and the sheets share their footprint.
      double fsrCourseWidth = 0.;

      // The last lead layer's front face and thickness. The left half
      // shares this plane, so it is recorded as the loop runs rather
      // than worked back out of fsrZ afterwards, which by then has
      // walked past the copper lining as well.
      double fsrLastLeadZ = 0.;
      double fsrLastLeadT = 0.;

      for (int L = 1; L <= fsrLayerN; ++L) {
        std::ostringstream key;
        key << "stm.FrontShieldingRight.layer" << L;
        const std::string kind = _config.getString(key.str());

        if (kind == "Pb") {
          ++fsrLeadLayer;

          // How deep the layer is: the bricks' own thickness, rotated. Read
          // from the first brick of the first course, since a layer is
          // one brick deep and they all lie the same way up.
          double thick = 0.;

          std::vector<BrickWallBrick> bricks;
          for (size_t c = 0; c < fsrCourse.size(); ++c) {
            std::ostringstream ckey;
            ckey << "stm.FrontShieldingRight.leadLayer" << fsrLeadLayer
                 << "Course" << c+1;
            std::vector<int> types;
            _config.getVectorInt(ckey.str(), types);
            if (types.empty()) {
              throw cet::exception("GEOM")
                << "STMMaker: " << ckey.str() << " is an empty course.\n";
            }

            // Which bores pass through the brick at a given position of
            // this course, read off the bore map. The holes are pinned
            // to the beam while brickEndX pins the courses, so the two
            // can be moved apart; layCourse checks each one lands on
            // the brick the map names.
            CourseBores fsrCourseBores;
            fsrCourseBores.bores           = &fsrBores;
            fsrCourseBores.wearAlongCourse = fsrWearAlong(fsrCourseDir);
            fsrCourseBores.wearAlongPitch  = fsrWearAlong(fsrPitchDir);
            fsrCourseBores.context         = "stm.FrontShieldingRight.leadBoreMap";
            fsrCourseBores.at = [&, c](int position) {
              std::vector<int> ids;
              for (auto & e : fsrBoreMap) {
                if (e.layer == fsrLeadLayer && e.course == int(c)+1
                    && e.position == position) {
                  ids.push_back(e.bore);
                  e.used = true;
                }
              }
              return ids;
            };

            // Butt the types end to end, flush at brickEndX: the course
            // runs toward -x, so that end IS the corner to start from.
            //
            // Unlike the walls the course plane here is absolute rather
            // than accumulated -- it comes off the baseplate chain as a
            // signed course number -- so layCourse is given the
            // course's lower face and re-centers each brick itself.
            // Backing off by the brick's own span rather than by
            // coursePitch keeps the center where it was asked for even
            // if a course is ever laid with thinner bricks.
            //
            // These are final positions: this section carries no
            // whole-structure shift, being pinned by brickEndX,
            // boreToBase and wallToCradleGap.
            const double courseHeight =
              brickSpan(_leadBrickDims, types[0], fsrBrickRot, fsrPitchDir,
                        ckey.str());

            const CLHEP::Hep3Vector fsrCourseStart(
              fsrBrickEndX,
              fsrCourseOffsetY + fsrCourse[c]*fsrPitch - courseHeight/2,
              fsrZ);

            const std::vector<BrickWallBrick> course =
              layCourse(types, fsrBrickOrient, fsrCourseStart,
                        fsrCourseDir, fsrPitchDir, fsrDepthDir,
                        _leadBrickDims, ckey.str(), 0, &thick,
                        fsrCourseBores);
            bricks.insert(bricks.end(), course.begin(), course.end());

            // The sheets share the courses' footprint.
            fsrCourseWidth = 0.;
            for (int t : types) fsrCourseWidth += fsrTypeWidth(t);
          }

          _FrontShieldingRightLeadLayers.push_back(
            BrickWall(_FrontShieldingRightBuild,
                      CLHEP::Hep3Vector(fsrBrickEndX, 0., fsrZ),
                      fsrCourseDir, fsrPitchDir, fsrDepthDir,
                      bricks, fsrBores));
          fsrLastLeadZ = fsrZ;
          fsrLastLeadT = thick;
          fsrZ += thick;
          _FrontShieldingRightLastLeadBackZ = fsrZ;

        } else {
          // A sheet layer: two plates split at the beam, sharing the
          // courses' footprint. The lower one carries the bores.
          //
          // The pair's top edge is level with the top face of the lead
          // stack beside it. That top face belongs to the highest lead
          // course, which is centred at fsrCourseOffsetY + nTop*fsrPitch
          // and is one pitch tall, so it sits half a pitch above that
          // center. The sheets then hang down from there by their own
          // heights, dyUpper and dyLower.
          const int fsrTopCourse =
            *std::max_element(fsrCourse.begin(), fsrCourse.end());
          const double top = fsrCourseOffsetY + fsrTopCourse*fsrPitch + fsrPitch/2;
          const double bot = top - (fsrSheetLo + fsrSheetUp);

          FrontShieldingRightSheet lower;
          lower.material    = fsrSheetMat;
          lower.orientation = "000";
          lower.halfDim     = CLHEP::Hep3Vector(fsrCourseWidth/2, fsrSheetLo/2, fsrSheetT/2);
          lower.center      = CLHEP::Hep3Vector(fsrBrickEndX - fsrCourseWidth/2,
                                                bot + fsrSheetLo/2,
                                                fsrZ + fsrSheetT/2);
          for (int i = 0; i < fsrBoreN; ++i) {
            lower.bores.push_back(i+1);
            lower.boreRadius.push_back(fsrBores[i].radius);
          }
          _FrontShieldingRightSheets.push_back(lower);

          FrontShieldingRightSheet upper;
          upper.material    = fsrSheetMat;
          upper.orientation = "000";
          upper.halfDim     = CLHEP::Hep3Vector(fsrCourseWidth/2, fsrSheetUp/2, fsrSheetT/2);
          upper.center      = CLHEP::Hep3Vector(fsrBrickEndX - fsrCourseWidth/2,
                                                bot + fsrSheetLo + fsrSheetUp/2,
                                                fsrZ + fsrSheetT/2);
          _FrontShieldingRightSheets.push_back(upper);

          fsrZ += fsrSheetT;
        }
      }

      // Every bore-map entry should have matched a brick by now.
      for (size_t i = 0; i < fsrBoreMap.size(); ++i) {
        FSRBoreMapEntry const & e = fsrBoreMap[i];
        if (!e.used) {
          throw cet::exception("GEOM")
            << "STMMaker: stm.FrontShieldingRight.leadBoreMap" << i+1
            << " = {" << e.layer << ", " << e.course << ", " << e.position
            << ", " << e.bore << "} names lead layer " << e.layer
            << ", course " << e.course << ", position " << e.position
            << ", but no brick is there (" << fsrLeadLayer << " lead layers, "
            << fsrCourse.size() << " courses each).\n";
        }
      }

      // The blocks the wall stands on: they span the courses' width and
      // their undersides land on the baseplate.
      _FrontShieldingRightBlockMaterial =
        _config.getString("stm.FrontShieldingRight.block.material");
      const double fsrBlockDx = _config.getDouble("stm.FrontShieldingRight.block.dx");
      const double fsrBlockDy = _config.getDouble("stm.FrontShieldingRight.block.dy");
      const double fsrBlockDz = _config.getDouble("stm.FrontShieldingRight.block.dz");
      _FrontShieldingRightBlockHalfDim =
        CLHEP::Hep3Vector(fsrBlockDx/2, fsrBlockDy/2, fsrBlockDz/2);

      const int fsrBlockN = _config.getInt("stm.FrontShieldingRight.block.n");
      const double fsrBlockY = -_STM_SSCboreToBase + fsrBlockDy/2;
      for (int i = 0; i < fsrBlockN; ++i) {
        _FrontShieldingRightBlockCenter.push_back(
          CLHEP::Hep3Vector(fsrBrickEndX - fsrBlockDx/2 - i*fsrBlockDx,
                            fsrBlockY,
                            fsrFrontZ + fsrBlockDz/2));
      }

      // The two pipes, each on its bore, starting at the wall front.
      const std::string fsrPipeMat = _config.getString("stm.FrontShieldingRight.pipe.material");
      const double fsrPipeRInLaBr  = _config.getDouble("stm.FrontShieldingRight.pipe.rInLaBr");
      const double fsrPipeRInHPGe  = _config.getDouble("stm.FrontShieldingRight.pipe.rInHPGe");
      const double fsrPipeROutLaBr = _config.getDouble("stm.FrontShieldingRight.pipe.rOutLaBr");
      const double fsrPipeROutHPGe = _config.getDouble("stm.FrontShieldingRight.pipe.rOutHPGe");
      const double fsrPipeLenLaBr  = _config.getDouble("stm.FrontShieldingRight.pipe.lenLaBr");
      const double fsrPipeLenHPGe  = _config.getDouble("stm.FrontShieldingRight.pipe.lenHPGe");
      for (int i = 0; i < fsrBoreN; ++i) {
        const bool laBr = (fsrBores[i].axis == "LaBr");
        FrontShieldingRightPipe p;
        p.material   = fsrPipeMat;
        p.rIn        = laBr ? fsrPipeRInLaBr  : fsrPipeRInHPGe;
        p.rOut       = laBr ? fsrPipeROutLaBr : fsrPipeROutHPGe;
        p.halfLength = (laBr ? fsrPipeLenLaBr : fsrPipeLenHPGe)/2;
        p.center     = CLHEP::Hep3Vector(fsrBores[i].center.x(),
                                         fsrBores[i].center.y(),
                                         fsrFrontZ + p.halfLength);
        _FrontShieldingRightPipes.push_back(p);
      }

      // The copper plate, placed by the anchor corner of its outline
      // rather than by a center.
      _FrontShieldingRightPlate.material =
        _config.getString("stm.FrontShieldingRight.copperLining.material");
      _config.getVectorDouble("stm.FrontShieldingRight.copperLining.UVerts",
                              _FrontShieldingRightPlate.uVerts);
      _config.getVectorDouble("stm.FrontShieldingRight.copperLining.VVerts",
                              _FrontShieldingRightPlate.vVerts);
      _FrontShieldingRightPlate.length =
        _config.getDouble("stm.FrontShieldingRight.copperLining.length");
      _FrontShieldingRightPlate.orientation =
        _config.getString("stm.FrontShieldingRight.copperLining.orientation");
      const double fsrCuLiningFromEnd =
        _config.getDouble("stm.FrontShieldingRight.copperLining.originFromBrickEnd");
      // The lining is not centered on the beam: its underside sits a
      // stated height above the baseplate, on the same chain that
      // fixes the bore height. The sweep runs along y and an extruded
      // solid is centered on its placement point, so the anchor's y is
      // the mid-plane -- bottom plus half the sweep -- not the edge.
      const double fsrCuLiningBottomToBase =
        _config.getDouble("stm.FrontShieldingRight.copperLining.bottomToBase");
      const double fsrCuLiningY =
        -_STM_SSCboreToBase + fsrCuLiningBottomToBase
        + _FrontShieldingRightPlate.length/2;
      _FrontShieldingRightPlate.anchor =
        CLHEP::Hep3Vector(fsrBrickEndX - fsrCuLiningFromEnd, fsrCuLiningY, fsrZ);
      // The plate's holes are the same bores as everything else on the
      // beamline, so it names them by id: whatever offset has moved a
      // bore off the collimator axis has already been applied there.
      const int fsrCuLiningHoles = _config.getInt("stm.FrontShieldingRight.copperLining.nHoles");
      if (fsrCuLiningHoles > fsrBoreN) {
        throw cet::exception("GEOM")
          << "STMMaker: stm.FrontShieldingRight.copperLining.nHoles is "
          << fsrCuLiningHoles << " but only " << fsrBoreN << " bores are defined.\n";
      }
      for (int i = 1; i <= fsrCuLiningHoles; ++i) {
        std::ostringstream key;
        key << "stm.FrontShieldingRight.copperLining.holeRadius" << i;
        _FrontShieldingRightPlate.bores.push_back(i);
        _FrontShieldingRightPlate.holeRadius.push_back(_config.getDouble(key.str()));
      }

      // The section's downstream face: the lining's front, plus however
      // much of the lining lies along the depth direction.
      {
        CLHEP::HepRotation liningRot(CLHEP::HepRotation::IDENTITY);
        OrientationResolver OR;
        OR.getRotationFromOrientation(liningRot, _FrontShieldingRightPlate.orientation);

        auto range = [](std::vector<double> const & v) {
          return v.empty() ? 0. : *std::max_element(v.begin(), v.end())
                                - *std::min_element(v.begin(), v.end());
        };
        // In the outline's own frame u and v span the polygon and the
        // sweep runs along its z.
        const CLHEP::Hep3Vector liningDim(range(_FrontShieldingRightPlate.uVerts),
                                          range(_FrontShieldingRightPlate.vVerts),
                                          _FrontShieldingRightPlate.length);
        _FrontShieldingRightBackZ =
          fsrZ + spanAlong(liningDim, liningRot, fsrDepthDir);
      }

      // ---- Front shielding, left ----------------------------------
      //
      // The half beside the collimator. Nothing here is bored and no
      // poly interleaves with a lead layer, so it is three brick
      // groups, four sheets and a prism, each placed against the
      // right half rather than measured on its own.
      //
      // Every x is a center measured back from brickEndX and so is
      // negative; every y is a height above the baseplate. The
      // per-bore offsets do not enter: they move the beam, and these
      // pieces are anchored by brickEndX and the baseplate.

      _FrontShieldingLeftBuild = _config.getBool("stm.FrontShieldingLeft.build");

      const std::string fslPolyMat =
        _config.getString("stm.FrontShieldingLeft.BPmaterial");

      auto fslSheet = [&](std::string const & name,
                          double dx, double dy, double dz,
                          double z, double bottomToBase) {
        FrontShieldingLeftSheet s;
        s.material = fslPolyMat;
        s.halfDim  = CLHEP::Hep3Vector(dx/2, dy/2, dz/2);
        s.center   = CLHEP::Hep3Vector(
          fsrBrickEndX + _config.getDouble("stm.FrontShieldingLeft." + name + ".fromEndX"),
          -_STM_SSCboreToBase + bottomToBase + dy/2,
          z + dz/2);
        return s;
      };

      // The side slab stands on the baseplate, spanning the depth of
      // the right half's first four layers from its front face.
      const double fslSideDx = _config.getDouble("stm.FrontShieldingLeft.sidePoly.dx");
      const double fslSideDy = _config.getDouble("stm.FrontShieldingLeft.sidePoly.dy");
      const double fslSideDz = _config.getDouble("stm.FrontShieldingLeft.sidePoly.dz");
      _FrontShieldingLeftSheets.push_back(
        fslSheet("sidePoly", fslSideDx, fslSideDy, fslSideDz, fsrFrontZ, 0.));

      // The outer sheet is the layer in front of the brick block, and
      // the block shares the plane of the right half's last lead
      // layer, so the sheet sits one thickness ahead of that plane.
      const double fslOuterDx = _config.getDouble("stm.FrontShieldingLeft.outerPoly.dx");
      const double fslOuterDy = _config.getDouble("stm.FrontShieldingLeft.outerPoly.dy");
      const double fslOuterDz = _config.getDouble("stm.FrontShieldingLeft.outerPoly.dz");
      _FrontShieldingLeftSheets.push_back(
        fslSheet("outerPoly", fslOuterDx, fslOuterDy, fslOuterDz,
                 fsrLastLeadZ - fslOuterDz, 0.));

      // The inner sheet and the edge sheet share a depth -- the layer
      // behind the last lead layer -- but not a height.
      const double fslInnerDx = _config.getDouble("stm.FrontShieldingLeft.innerPoly.dx");
      const double fslInnerDy = _config.getDouble("stm.FrontShieldingLeft.innerPoly.dy");
      const double fslInnerDz = _config.getDouble("stm.FrontShieldingLeft.innerPoly.dz");
      const double fslInnerBase =
        _config.getDouble("stm.FrontShieldingLeft.innerPoly.bottomToBase");
      const double fslInnerZ = fsrLastLeadZ + fsrLastLeadT;
      _FrontShieldingLeftSheets.push_back(
        fslSheet("innerPoly", fslInnerDx, fslInnerDy, fslInnerDz,
                 fslInnerZ, fslInnerBase));

      const double fslEdgeDx = _config.getDouble("stm.FrontShieldingLeft.edgePoly.dx");
      const double fslEdgeDy = _config.getDouble("stm.FrontShieldingLeft.edgePoly.dy");
      const double fslEdgeDz = _config.getDouble("stm.FrontShieldingLeft.edgePoly.dz");
      const double fslEdgeBase =
        _config.getDouble("stm.FrontShieldingLeft.edgePoly.bottomToBase");
      _FrontShieldingLeftSheets.push_back(
        fslSheet("edgePoly", fslEdgeDx, fslEdgeDy, fslEdgeDz,
                 fslInnerZ, fslEdgeBase));

      // A brick group: one BrickWall, its courses butted end to end
      // from a corner. The same expansion serves all three groups
      // here, which differ only in where they start and which way
      // their courses run. Orientation is per course, since these
      // groups mix bricks lying flat with bricks stood on end.
      auto fslBrickGroup = [&](std::vector<std::vector<int> > const & courses,
                               std::vector<std::string> const & orientations,
                               CLHEP::Hep3Vector const & origin,
                               CLHEP::Hep3Vector const & courseDir,
                               CLHEP::Hep3Vector const & pitchDir,
                               CLHEP::Hep3Vector const & depthDir,
                               std::string const & what) {
        std::vector<BrickWallBrick> bricks;
        double pitchAt = 0.;
        for (size_t c = 0; c < courses.size(); ++c) {
          std::ostringstream ckey;
          ckey << "FrontShieldingLeft " << what << " course " << c+1;

          // The step is read per course rather than once, because
          // these groups mix brick sizes from course to course.
          double step = 0.;
          const std::vector<BrickWallBrick> course =
            layCourse(courses[c], orientations[c],
                      origin + pitchDir*pitchAt,
                      courseDir, pitchDir, depthDir,
                      _leadBrickDims, ckey.str(), &step);
          bricks.insert(bricks.end(), course.begin(), course.end());
          pitchAt += step;
        }
        return BrickWall(_FrontShieldingLeftBuild, origin,
                         courseDir, pitchDir, depthDir,
                         bricks, std::vector<BrickWallBore>());
      };

      // The nine-brick block. Its columns run up y, so a "course"
      // here is a column: the course direction is +y and the pitch
      // steps along -x, away from the right half. It stands on the
      // baseplate.
      //
      // Its +x face is flush with the side slab's +x face, which is the
      // right half's -x course edge, but going through the slab is
      // what keeps them together if the slab is ever moved or
      // retuned. The slab stands alongside the block in z, spanning
      // the right half's first four layers.
      const double fslBlockEndX =
        fsrBrickEndX + _config.getDouble("stm.FrontShieldingLeft.sidePoly.fromEndX")
        + fslSideDx/2;
      {
        std::vector<std::vector<int> > columns;
        std::vector<std::string>       orientations;
        const int n = _config.getInt("stm.FrontShieldingLeft.blockColumnN");
        for (int i = 1; i <= n; ++i) {
          std::ostringstream key, okey;
          key  << "stm.FrontShieldingLeft.blockColumn" << i;
          okey << "stm.FrontShieldingLeft.blockColumn" << i << "Orientation";
          std::vector<int> types;
          _config.getVectorInt(key.str(), types);
          columns.push_back(types);
          orientations.push_back(_config.getString(okey.str()));
        }
        _FrontShieldingLeftBrickGroups.push_back(
          fslBrickGroup(columns, orientations,
                        CLHEP::Hep3Vector(fslBlockEndX, -_STM_SSCboreToBase,
                                          fsrLastLeadZ),
                        CLHEP::Hep3Vector(0., 1., 0.),   // a column runs up
                        CLHEP::Hep3Vector(-1., 0., 0.),  // columns step -x
                        CLHEP::Hep3Vector(0., 0., 1.),
                        "outer block"));
      }

      // The two bricks beside the edge sheet. They are level with it
      // rather than stacked on it, so they take its height, and they
      // sit one layer along +z from it.
      {
        std::vector<int> types;
        _config.getVectorInt("stm.FrontShieldingLeft.edgeCourse", types);
        const std::string orient =
          _config.getString("stm.FrontShieldingLeft.edgeCourseOrientation");
        // Flush with the sheet's +x end, running -x, the same way the
        // block and the right half's courses run.
        const double xStart =
          fsrBrickEndX + _config.getDouble("stm.FrontShieldingLeft.edgePoly.fromEndX")
          + fslEdgeDx/2;
        _FrontShieldingLeftBrickGroups.push_back(
          fslBrickGroup(std::vector<std::vector<int> >(1, types),
                        std::vector<std::string>(1, orient),
                        CLHEP::Hep3Vector(xStart,
                                          -_STM_SSCboreToBase + fslEdgeBase,
                                          fslInnerZ + fslEdgeDz),
                        CLHEP::Hep3Vector(-1., 0., 0.),
                        CLHEP::Hep3Vector(0., 1., 0.),
                        CLHEP::Hep3Vector(0., 0., 1.),   // +z of the sheet
                        "edge course"));
      }

      // The grid behind the inner sheet: courses running -x from the
      // sheet's +x end, stacking up from the baseplate.
      {
        std::vector<std::vector<int> > courses;
        std::vector<std::string>       orientations;
        const int n = _config.getInt("stm.FrontShieldingLeft.innerGridCourseN");
        const std::string orient =
          _config.getString("stm.FrontShieldingLeft.innerGridOrientation");
        for (int i = 1; i <= n; ++i) {
          std::ostringstream key;
          key << "stm.FrontShieldingLeft.innerGridCourse" << i;
          std::vector<int> types;
          _config.getVectorInt(key.str(), types);
          courses.push_back(types);
          orientations.push_back(orient);
        }
        const double xStart =
          fsrBrickEndX + _config.getDouble("stm.FrontShieldingLeft.innerPoly.fromEndX")
          + fslInnerDx/2;
        _FrontShieldingLeftBrickGroups.push_back(
          fslBrickGroup(courses, orientations,
                        CLHEP::Hep3Vector(xStart,
                                          -_STM_SSCboreToBase + fslInnerBase,
                                          fslInnerZ + fslInnerDz),
                        CLHEP::Hep3Vector(-1., 0., 0.),
                        CLHEP::Hep3Vector(0., 1., 0.),
                        CLHEP::Hep3Vector(0., 0., 1.),   // +z of the sheet
                        "inner grid"));
      }

      // The prism, placed by its right angle rather than by a center.
      _FrontShieldingLeftPrism.material =
        _config.getString("stm.FrontShieldingLeft.triangle.material");
      _config.getVectorDouble("stm.FrontShieldingLeft.triangle.UVerts",
                              _FrontShieldingLeftPrism.uVerts);
      _config.getVectorDouble("stm.FrontShieldingLeft.triangle.VVerts",
                              _FrontShieldingLeftPrism.vVerts);
      _FrontShieldingLeftPrism.length =
        _config.getDouble("stm.FrontShieldingLeft.triangle.length");
      _FrontShieldingLeftPrism.orientation =
        _config.getString("stm.FrontShieldingLeft.triangle.orientation");
      // The right angle butts the inner sheet's +x and -z faces, and
      // the prism rests on the floor, its bottom bottomToBase above the
      // baseplate. An extruded solid is centered on its placement point,
      // so the anchor's y is half a sweep above that bottom.
      {
        const double innerX =
          fsrBrickEndX + _config.getDouble("stm.FrontShieldingLeft.innerPoly.fromEndX");
        const double prismBottom =
          -_STM_SSCboreToBase
          + _config.getDouble("stm.FrontShieldingLeft.triangle.bottomToBase");
        _FrontShieldingLeftPrism.anchor =
          CLHEP::Hep3Vector(innerX + fslInnerDx/2,
                            prismBottom + _FrontShieldingLeftPrism.length/2,
                            fslInnerZ);
      }
    }
    else{
      _FrontShieldingBuild   = _config.getBool(  "stm.FrontShielding.build");
      _FrontSHeightofRoom    = _config.getDouble("stm.FrontShielding.HeightofRoom");
      _FrontStungstenlength  = _config.getDouble("stm.FrontShielding.tungstenlength");
      _FrontStungstendepth   = _config.getDouble("stm.FrontShielding.tungstendepth");
      _FrontSleaddepth1      = _config.getDouble("stm.FrontShielding.leaddepth1");
      _FrontSleaddepth2      = _config.getDouble("stm.FrontShielding.leaddepth2");
      _FrontSaluminumdepth   = _config.getDouble("stm.FrontShielding.aluminumdepth");
      _FrontScopperdepth     = _config.getDouble("stm.FrontShielding.copperdepth");
      _FrontSBPdepth         = _config.getDouble("stm.FrontShielding.BPdepth");
      _FrontSfPb_lengthL     = _config.getDouble("stm.FrontShielding.fPb_lengthL");
      _FrontSfPb_lengthR     = _config.getDouble("stm.FrontShielding.fPb_lengthR");
      _FrontSGapForTop       = _config.getDouble("stm.FrontShielding.GapForTop");
      _FrontSLeakForSSC      = _config.getDouble("stm.FrontShielding.LeakForSSC");
      _FrontSCopperL         = _config.getDouble("stm.FrontShielding.CopperL");
      _FrontS_H              = _config.getDouble("stm.FrontShielding.FrontS_H");
      _FrontSHole_r          = _config.getDouble("stm.FrontShielding.FrontSHole_r");
    }


    _HPGeBuild                = _config.getBool("stm.HPGe.build");
    _HPGecrystalMaterial      = _config.getString("stm.HPGe.crystalMaterial");
    _HPGeholeMaterial         = _config.getString("stm.HPGe.holeMaterial");
    _HPGewindowMaterial       = _config.getString("stm.HPGe.windowMaterial");
    _HPGewallMaterial         = _config.getString("stm.HPGe.wallMaterial");
    _HPGecapsuleMaterial      = _config.getString("stm.HPGe.capsuleMaterial");
    _HPGeEndcapR              = _config.getDouble("stm.HPGe.EndcapR");
    _HPGeEndcapL              = _config.getDouble("stm.HPGe.EndcapL");
    _HPGeCrystalR             = _config.getDouble("stm.HPGe.CrystalR");
    _HPGeCrystalL             = _config.getDouble("stm.HPGe.CrystalL");
    _HPGeZ_HPGe               = _config.getDouble("stm.HPGe.Z_HPGe");
    _HPGeHoleR                = _config.getDouble("stm.HPGe.HoleR");
    _HPGeHoleL                = _config.getDouble("stm.HPGe.HoleL");
    _HPGeCapsule_Wallthick    = _config.getDouble("stm.HPGe.Capsule_Wallthick");
    _HPGeCapsule_Windowthick  = _config.getDouble("stm.HPGe.Capsule_Windowthick");
    _HPGeCapsule_Endthick     = _config.getDouble("stm.HPGe.Capsule_Endthick");
    _HPGeCapsule_Walllength   = _config.getDouble("stm.HPGe.Capsule_Walllength");
    _HPGeWindowD              = _config.getDouble("stm.HPGe.WindowD");
    _HPGeEndcapD              = _config.getDouble("stm.HPGe.EndcapD");
    _HPGeAirD                 = _config.getDouble("stm.HPGe.AirD");
    _HPGeoffset_HPGe          = _config.getDouble("stm.HPGe.offset_HPGe");

    _LaBrBuild                = _config.getBool("stm.LaBr.build");
    _LaBrcrystalMaterial      = _config.getString("stm.LaBr.crystalMaterial");
    _LaBrwindowMaterial       = _config.getString("stm.LaBr.windowMaterial");
    _LaBrwallMaterial         = _config.getString("stm.LaBr.wallMaterial");
    _LaBrEndcapR              = _config.getDouble("stm.LaBr.EndcapR");
    _LaBrEndcapL              = _config.getDouble("stm.LaBr.EndcapL");
    _LaBrCrystalR             = _config.getDouble("stm.LaBr.CrystalR");
    _LaBrCrystalL             = _config.getDouble("stm.LaBr.CrystalL");
    _LaBrZ_LaBr               = _config.getDouble("stm.LaBr.Z_LaBr");
    _LaBrWindowD              = _config.getDouble("stm.LaBr.WindowD");
    _LaBrEndcapD              = _config.getDouble("stm.LaBr.EndcapD");
    _LaBrAirD                 = _config.getDouble("stm.LaBr.AirD");
    _LaBroffset_LaBr          = _config.getDouble("stm.LaBr.offset_LaBr");


    // The shield house's reference corner, where the walls meet the
    // front shielding: the +x end of its courses, and the plane a
    // wall's own offsets measure from. Resolved here, before any wall
    // is expanded, so the walls key off one point instead of each
    // rebuilding it -- and so none of them has to read another wall's
    // keys to find it.
    if (_handstacked) {
      _ShieldHouseRefX = _FrontShieldingRightBrickEndX;
      _ShieldHouseRefZ = _FrontShieldingRightFrontZ
                       + _config.getDouble("stm.RightWall.refFromFrontZ");
    }

    // The updated bottom wall replaces BottomShielding. The two read
    // different keys and a geometry file carries only one set.
    if (_handstacked) {
      // ---- Bottom wall ------------------------------------------------
      //
      // The floor, including the baseplate the whole house stands on.
      // Five layers worked UPWARD in +y above the plate: a poly prism,
      // a lead layer, the same prism again, a second lead layer, and
      // two copper prisms.
      //
      // The reference is 5 in inboard of the house reference in x, at
      // that reference's z, and at the TOP of the baseplate in y --
      // the datum every wall stands on.

      _BottomWallBuild = _config.getBool("stm.BottomWall.build");

      const double bwRefX =
        _ShieldHouseRefX + _config.getDouble("stm.BottomWall.refFromHouseX");
      const double bwRefZ = _ShieldHouseRefZ;
      const double bwRefY = -_STM_SSCboreToBase;

      const std::string bwPolyMat = _config.getString("stm.BottomWall.BPmaterial");

      // ---- the baseplate ----------------------------------------------
      {
        const double dx = _config.getDouble("stm.BottomWall.basePlate.dx");
        const double dy = _config.getDouble("stm.BottomWall.basePlate.dy");
        const double dz = _config.getDouble("stm.BottomWall.basePlate.dz");

        _BottomWallPlate.material = _config.getString("stm.BottomWall.basePlate.material");
        _BottomWallPlate.halfDim  = CLHEP::Hep3Vector(dx/2, dy/2, dz/2);
        // The reference is the plate's TOP, so the plate hangs below it
        // and its center is half a thickness down.
        _BottomWallPlate.center   = CLHEP::Hep3Vector(
          bwRefX + _config.getDouble("stm.BottomWall.basePlate.centerFromRefX"),
          bwRefY - dy/2,
          bwRefZ + _config.getDouble("stm.BottomWall.basePlate.centerFromRefZ"));
      }

      // Courses pitch along -x from the reference and run along +z.
      const CLHEP::Hep3Vector bwCourseDir( 0., 0., 1.);
      const CLHEP::Hep3Vector bwPitchDir (-1., 0., 0.);
      const CLHEP::Hep3Vector bwDepthDir ( 0., 1., 0.);  // layers stack up

      // Walk the layers upward from the top of the baseplate. Each
      // starts at the face of the one before it. As in the side walls
      // the cursor is a point advanced along depthDir, so the
      // placements below read the same way in all four walls whichever
      // way a given one stacks.
      const int bwLayerN = _config.getInt("stm.BottomWall.layerN");
      CLHEP::Hep3Vector bwAt(bwRefX, bwRefY, bwRefZ);
      int bwLeadLayer = 0;
      int bwPolyLayer = 0;

      for (int L = 1; L <= bwLayerN; ++L) {
        std::ostringstream key;
        key << "stm.BottomWall.layer" << L;
        const std::string kind = _config.getString(key.str());

        if (kind == "BP") {
          // The poly L, the same shape both times it appears.
          ++bwPolyLayer;

          const double t = _config.getDouble("stm.BottomWall.polyPrism.thickness");

          std::ostringstream name;
          name << "BPBottomWallPrism" << bwPolyLayer << "PV";

          BottomWallPrism p;
          p.name     = name.str();
          p.material = bwPolyMat;
          _config.getVectorDouble("stm.BottomWall.polyPrism.UVerts", p.uVerts);
          _config.getVectorDouble("stm.BottomWall.polyPrism.VVerts", p.vVerts);
          p.length      = t;
          p.orientation = _config.getString("stm.BottomWall.polyPrism.orientation");
          // The cap's origin is at the reference in x and z; the sweep
          // is centered on the placement point, so the anchor sits half
          // a thickness into this layer.
          p.anchor = bwAt + bwDepthDir*(t/2);
          _BottomWallPrisms.push_back(p);

          bwAt += bwDepthDir*t;

        } else if (kind == "Pb") {
          ++bwLeadLayer;

          // Each course states its own orientation, since they are not
          // all the same width, so no layer-wide one is passed.
          std::ostringstream pkey;
          pkey << "stm.BottomWall.leadLayer" << bwLeadLayer;

          double thick = 0.;
          std::vector<BrickWallBrick> bricks =
            layLayer(_config, pkey.str(), bwAt,
                     bwCourseDir, bwPitchDir, bwDepthDir,
                     _leadBrickDims, "", &thick);

          // Then the few that break the course pattern. This is the
          // only wall that has any.
          const std::vector<BrickWallBrick> strays =
            layStrays(_config, pkey.str(), bwAt, bwDepthDir,
                      _leadBrickDims, thick);
          bricks.insert(bricks.end(), strays.begin(), strays.end());

          _BottomWallLeadLayers.push_back(
            BrickWall(_BottomWallBuild, bwAt,
                      bwCourseDir, bwPitchDir, bwDepthDir,
                      bricks, std::vector<BrickWallBore>()));

          bwAt += bwDepthDir*thick;

        } else if (kind == "Cu") {
          // Two prisms, each with its own outline and its own offsets.
          const std::string mat =
            _config.getString("stm.BottomWall.copperPrism.material");
          const double t =
            _config.getDouble("stm.BottomWall.copperPrism.thickness");
          const std::string orient =
            _config.getString("stm.BottomWall.copperPrism.orientation");
          // These do not sit on the layer below: their height above the
          // baseplate is stated, and the anchor is the mid-plane of the
          // sweep rather than a face.
          const double toBase =
            _config.getDouble("stm.BottomWall.copperPrism.anchorToBase");

          for (int i = 1; i <= 2; ++i) {
            std::ostringstream base, name;
            base << "stm.BottomWall.copperPrism" << i;
            name << "CopperBottomWallPrism" << i << "PV";

            BottomWallPrism p;
            p.name     = name.str();
            p.material = mat;
            _config.getVectorDouble(base.str() + ".UVerts", p.uVerts);
            _config.getVectorDouble(base.str() + ".VVerts", p.vVerts);
            p.length      = t;
            p.orientation = orient;
            // Measured from the reference, not from the cursor, per
            // anchorToBase above.
            p.anchor = CLHEP::Hep3Vector(
              bwRefX + _config.getDouble(base.str() + ".fromRefX"),
              bwRefY + toBase,
              bwRefZ + _config.getDouble(base.str() + ".fromRefZ"));
            _BottomWallPrisms.push_back(p);
          }

          bwAt += bwDepthDir*t;

        } else {
          throw cet::exception("GEOM")
            << "STMMaker: " << key.str() << " is \"" << kind
            << "\", which is not one of BP, Pb or Cu.\n";
        }
      }

    } else {
      _BottomShieldingBuild  = _config.getBool("stm.BottomShielding.build");
      _BottomSfloor_Zlength  = _config.getDouble("stm.BottomShielding.floor_Zlength");
      _BottomSFront_LB       = _config.getDouble("stm.BottomShielding.Front_LB");
      _BottomSFront_LB_inner = _config.getDouble("stm.BottomShielding.Front_LB_inner");
      _BottomSleaddepth      = _config.getDouble("stm.BottomShielding.leaddepth");
      _BottomScopperdepth    = _config.getDouble("stm.BottomShielding.copperdepth");
      _BottomSBPdepth        = _config.getDouble("stm.BottomShielding.BPdepth");
    }


    // The updated left wall replaces LeftShielding, as above.
    if (_handstacked) {
      // ---- Left wall --------------------------------------------------
      //
      // The mirror of the right wall: the same five kinds of layer and
      // the same sheet shapes, without the L-shaped edge piece.
      //
      // Layers are listed outside to inside and walk +x from the
      // reference, which is the -x face of the outermost layer, so every
      // piece lies on the same side of it.
      //
      // In z the reference is the -z end, shared with the right wall.
      // The sheets and both lead layers run +z from it; the copper is
      // the exception, sitting flush at the far +z end instead.

      _LeftWallBuild = _config.getBool("stm.LeftWall.build");

      const double lwRefX =
        _ShieldHouseRefX + _config.getDouble("stm.LeftWall.refFromBrickEndX");
      // The walls all start together in z, at the house reference.
      const double lwRefZ  = _ShieldHouseRefZ;
      const double lwBaseY = -_STM_SSCboreToBase;

      const std::string lwPolyMat = _config.getString("stm.LeftWall.BPmaterial");

      // Courses run +z from the reference and stack up from the
      // baseplate; the layers themselves build inward along +x.
      const CLHEP::Hep3Vector lwCourseDir(0., 0., 1.);
      const CLHEP::Hep3Vector lwPitchDir (0., 1., 0.);
      const CLHEP::Hep3Vector lwDepthDir (1., 0., 0.);

      const std::string lwBrickOrient =
        _config.getString("stm.LeftWall.brickOrientation");

      // How far the wall runs in z. The sheets and the lead layers all
      // share it, and the copper is placed flush against its far end,
      // so it is read once from the poly rather than per piece.
      const double lwRunZ = _config.getDouble("stm.LeftWall.longwallPoly.dz");

      // Walk the layers inward, each starting at the face of the one
      // before it. As in the right wall the cursor is a point advanced
      // along depthDir, so none of the placements below carries the
      // sign of the walk even though this wall stacks the other way.
      const int lwLayerN = _config.getInt("stm.LeftWall.layerN");
      CLHEP::Hep3Vector lwAt(lwRefX, lwBaseY, lwRefZ);
      int lwLeadLayer = 0;
      int lwPolyPair = 0;

      for (int L = 1; L <= lwLayerN; ++L) {
        std::ostringstream key;
        key << "stm.LeftWall.layer" << L;
        const std::string kind = _config.getString(key.str());

        if (kind == "Pb") {
          ++lwLeadLayer;

          // The layer is one brick deep, and the courses stack by the
          // brick's own span -- both follow from the rotation. One
          // orientation for the whole layer, as in the right wall.
          std::ostringstream pkey;
          pkey << "stm.LeftWall.leadLayer" << lwLeadLayer;

          double thick = 0.;
          const std::vector<BrickWallBrick> bricks =
            layLayer(_config, pkey.str(), lwAt,
                     lwCourseDir, lwPitchDir, lwDepthDir,
                     _leadBrickDims, lwBrickOrient, &thick);

          _LeftWallLeadLayers.push_back(
            BrickWall(_LeftWallBuild, lwAt,
                      lwCourseDir, lwPitchDir, lwDepthDir,
                      bricks, std::vector<BrickWallBore>()));
          lwAt += lwDepthDir*thick;

        } else if (kind == "BP") {
          // A longwall standing on the baseplate and a top edge above
          // it, meeting with no gap.
          ++lwPolyPair;

          const double lgDx = _config.getDouble("stm.LeftWall.longwallPoly.dx");
          const double lgDy = _config.getDouble("stm.LeftWall.longwallPoly.dy");
          const double lgDz = _config.getDouble("stm.LeftWall.longwallPoly.dz");
          const double lgToBase =
            _config.getDouble("stm.LeftWall.longwallPoly.bottomToBase");

          std::ostringstream lgName;
          lgName << "BPLeftWallLongwall" << lwPolyPair << "PV";

          WallSheet longwall;
          longwall.name     = lgName.str();
          longwall.material = lwPolyMat;
          longwall.halfDim  = CLHEP::Hep3Vector(lgDx/2, lgDy/2, lgDz/2);
          longwall.center   = lwAt + lwDepthDir*(lgDx/2)
                            + CLHEP::Hep3Vector(0., lgToBase + lgDy/2, lgDz/2);
          _LeftWallSheets.push_back(longwall);

          const double teDx = _config.getDouble("stm.LeftWall.topEdgePoly.dx");
          const double teDy = _config.getDouble("stm.LeftWall.topEdgePoly.dy");
          const double teDz = _config.getDouble("stm.LeftWall.topEdgePoly.dz");
          const double teToBase =
            _config.getDouble("stm.LeftWall.topEdgePoly.bottomToBase");

          std::ostringstream teName;
          teName << "BPLeftWallTopEdge" << lwPolyPair << "PV";

          WallSheet top;
          top.name     = teName.str();
          top.material = lwPolyMat;
          top.halfDim  = CLHEP::Hep3Vector(teDx/2, teDy/2, teDz/2);
          top.center   = lwAt + lwDepthDir*(teDx/2)
                       + CLHEP::Hep3Vector(0., teToBase + teDy/2, teDz/2);
          _LeftWallSheets.push_back(top);

          // Both are one plate thick, so either gives the step.
          lwAt += lwDepthDir*lgDx;

        } else if (kind == "Cu") {
          const double dx = _config.getDouble("stm.LeftWall.copperSheet.dx");
          const double dy = _config.getDouble("stm.LeftWall.copperSheet.dy");
          const double dz = _config.getDouble("stm.LeftWall.copperSheet.dz");
          const double toBase =
            _config.getDouble("stm.LeftWall.copperSheet.bottomToBase");

          // The one piece that does not start at the reference: it runs
          // the wall's length downstream and sits flush at the far +z
          // end, so its center is half its own depth back from there.
          WallSheet s;
          s.name     = "CopperLeftWallSheetPV";
          s.material = _config.getString("stm.LeftWall.copperSheet.material");
          s.halfDim  = CLHEP::Hep3Vector(dx/2, dy/2, dz/2);
          s.center   = lwAt + lwDepthDir*(dx/2)
                     + CLHEP::Hep3Vector(0., toBase + dy/2, lwRunZ - dz/2);
          _LeftWallSheets.push_back(s);

          lwAt += lwDepthDir*dx;

        } else {
          throw cet::exception("GEOM")
            << "STMMaker: " << key.str() << " is \"" << kind
            << "\", which is not one of Pb, BP or Cu.\n";
        }
      }

    } else {
      _LeftShieldingBuild   = _config.getBool("stm.LeftShielding.build");
      _LeftS_Length         = _config.getDouble("stm.LeftShielding.Length");
      _LeftSleaddepth       = _config.getDouble("stm.LeftShielding.leaddepth");
      _LeftScopperdepth     = _config.getDouble("stm.LeftShielding.copperdepth");
      _LeftSBPdepth         = _config.getDouble("stm.LeftShielding.BPdepth");
      _LeftSXmin            = _config.getDouble("stm.LeftShielding.Left_Xmin");
    }

    // The updated right wall replaces RightShielding, as above.
    if (_handstacked) {
      // ---- Right wall -------------------------------------------------
      //
      // The wall along +x, running downstream from the front shielding.
      // Five layers worked inward along -x from one reference vertex,
      // which is
      //
      //   - refFromFrontZ downstream of the front shielding's front
      //     face, flush with the +z side of its last lead layer,
      //   - at that layer's +x, -y corner.
      //
      // Courses run along +z and the layers stack along -x.

      _RightWallBuild = _config.getBool("stm.RightWall.build");

      const double rwRefFromFront =
        _config.getDouble("stm.RightWall.refFromFrontZ");
      // The house reference is built from this same key, so the two
      // agree by construction; it is named here for the messages below.
      const double rwRefZ = _ShieldHouseRefZ;

      // The reference is meant to land on the back of the front
      // shielding's last lead layer, a plane the layer stack also
      // gives. The wall runs downstream from the reference, so a
      // reference upstream of that plane drives the wall INTO the
      // front shielding; downstream of it only opens a gap.
      {
        const double dimTolerance = 0.001;  // mm
        const double slack = rwRefZ - _FrontShieldingRightLastLeadBackZ;
        if (slack < -dimTolerance) {
          throw cet::exception("GEOM")
            << "STMMaker: stm.RightWall.refFromFrontZ puts the right wall's"
            << " reference at z = " << rwRefZ << " mm, "
            << -slack << " mm INSIDE the front shielding, whose last lead"
            << " layer ends at " << _FrontShieldingRightLastLeadBackZ
            << " mm.\n"
            << "The two would overlap. Check refFromFrontZ ("
            << rwRefFromFront << " mm) against the front shielding's"
            << " layer sequence.\n";
        }
        else if (slack > dimTolerance) {
          mf::LogWarning("GEOM")
            << "STMMaker: the right wall starts " << slack
            << " mm downstream of the front shielding's last lead layer"
            << " (wall reference z = " << rwRefZ << " mm, layer ends at "
            << _FrontShieldingRightLastLeadBackZ << " mm), leaving a gap."
            << " Intentional if refFromFrontZ (" << rwRefFromFront
            << " mm) was set that way.\n";
        }
        else if (_verbosityLevel > 0) {
          mf::LogInfo("GEOM")
            << "STMMaker: the right wall butts the front shielding's last"
            << " lead layer (<" << dimTolerance << " mm), at z = " << rwRefZ
            << " mm.\n";
        }
      }

      // The wall's +x face: the front shielding's courses end there.
      const double rwRefX  = _ShieldHouseRefX;
      const double rwBaseY = -_STM_SSCboreToBase;

      const std::string rwPolyMat = _config.getString("stm.RightWall.BPmaterial");

      // Courses run +z from the reference and stack up from the
      // baseplate; the layers themselves step inward along -x.
      const CLHEP::Hep3Vector rwCourseDir( 0., 0., 1.);
      const CLHEP::Hep3Vector rwPitchDir ( 0., 1., 0.);
      const CLHEP::Hep3Vector rwDepthDir (-1., 0., 0.);

      const std::string rwBrickOrient =
        _config.getString("stm.RightWall.brickOrientation");

      // Walk the layers inward. Each starts at the face of the one
      // before it, so a layer's own thickness is the only thing that
      // separates it from the next.
      //
      // The cursor is a point rather than an x, and it advances along
      // depthDir. That way nothing here carries the sign of the walk:
      // a piece is placed half its own thickness along depthDir from
      // the cursor and the cursor moves on by the whole thickness,
      // reading the same in this wall as in the left one that stacks
      // the other way.
      const int rwLayerN = _config.getInt("stm.RightWall.layerN");
      CLHEP::Hep3Vector rwAt(rwRefX, rwBaseY, rwRefZ);
      int rwLeadLayer = 0;

      for (int L = 1; L <= rwLayerN; ++L) {
        std::ostringstream key;
        key << "stm.RightWall.layer" << L;
        const std::string kind = _config.getString(key.str());

        if (kind == "Pb") {
          ++rwLeadLayer;

          // The layer is one brick deep, and the courses stack by the
          // brick's own span -- both follow from the rotation, so
          // neither is given in the config. One orientation for the
          // whole layer, so the courses state none of their own.
          std::ostringstream pkey;
          pkey << "stm.RightWall.leadLayer" << rwLeadLayer;

          double thick = 0.;
          const std::vector<BrickWallBrick> bricks =
            layLayer(_config, pkey.str(), rwAt,
                     rwCourseDir, rwPitchDir, rwDepthDir,
                     _leadBrickDims, rwBrickOrient, &thick);

          _RightWallLeadLayers.push_back(
            BrickWall(_RightWallBuild, rwAt,
                      rwCourseDir, rwPitchDir, rwDepthDir,
                      bricks, std::vector<BrickWallBore>()));
          rwAt += rwDepthDir*thick;

        } else if (kind == "BPouter") {
          // The L-shaped edge piece, placed by the corner where it
          // meets the front shielding, and the sheet above it.
          _RightWallEdgePrism.material = rwPolyMat;
          _config.getVectorDouble("stm.RightWall.edgePoly.UVerts",
                                  _RightWallEdgePrism.uVerts);
          _config.getVectorDouble("stm.RightWall.edgePoly.VVerts",
                                  _RightWallEdgePrism.vVerts);
          _RightWallEdgePrism.length =
            _config.getDouble("stm.RightWall.edgePoly.thickness");
          _RightWallEdgePrism.orientation =
            _config.getString("stm.RightWall.edgePoly.orientation");
          // The sweep is centered on the placement point, so the anchor
          // sits half a thickness inboard of this layer's outer face.
          // This layer is the one exception to the inward walk: it lies
          // on the far side of the cursor from the rest, hence the step
          // AGAINST depthDir.
          _RightWallEdgePrism.anchor =
            rwAt - rwDepthDir*(_RightWallEdgePrism.length/2);

          const double dx = _config.getDouble("stm.RightWall.outerTopPoly.dx");
          const double dy = _config.getDouble("stm.RightWall.outerTopPoly.dy");
          const double dz = _config.getDouble("stm.RightWall.outerTopPoly.dz");
          const double toBase =
            _config.getDouble("stm.RightWall.outerTopPoly.bottomToBase");

          // Its -z end is flush with the +z side of the L's upright,
          // which stands at the reference, so its z follows from that
          // rather than being given.
          WallSheet s;
          s.name     = "BPRightWallOuterTopPV";
          s.material = rwPolyMat;
          s.halfDim  = CLHEP::Hep3Vector(dx/2, dy/2, dz/2);
          s.center   = rwAt - rwDepthDir*(dx/2)
                     + CLHEP::Hep3Vector(0., toBase + dy/2, dz/2);
          _RightWallSheets.push_back(s);

          // The cursor is not advanced here: this layer sits on the
          // other side of the reference point from the inward walk.

        } else if (kind == "BPinner") {
          // The longwall below and the top sheet above it, sharing this
          // layer's plane and the reference as their -z end.
          const double lwDx = _config.getDouble("stm.RightWall.longwallPoly.dx");
          const double lwDy = _config.getDouble("stm.RightWall.longwallPoly.dy");
          const double lwDz = _config.getDouble("stm.RightWall.longwallPoly.dz");
          const double lwToBase =
            _config.getDouble("stm.RightWall.longwallPoly.bottomToBase");

          WallSheet longwall;
          longwall.name     = "BPRightWallLongwallPV";
          longwall.material = rwPolyMat;
          longwall.halfDim  = CLHEP::Hep3Vector(lwDx/2, lwDy/2, lwDz/2);
          longwall.center   = rwAt + rwDepthDir*(lwDx/2)
                            + CLHEP::Hep3Vector(0., lwToBase + lwDy/2, lwDz/2);
          _RightWallSheets.push_back(longwall);

          const double tpDx = _config.getDouble("stm.RightWall.innerTopPoly.dx");
          const double tpDy = _config.getDouble("stm.RightWall.innerTopPoly.dy");
          const double tpDz = _config.getDouble("stm.RightWall.innerTopPoly.dz");
          const double tpToBase =
            _config.getDouble("stm.RightWall.innerTopPoly.bottomToBase");

          WallSheet top;
          top.name     = "BPRightWallInnerTopPV";
          top.material = rwPolyMat;
          top.halfDim  = CLHEP::Hep3Vector(tpDx/2, tpDy/2, tpDz/2);
          top.center   = rwAt + rwDepthDir*(tpDx/2)
                       + CLHEP::Hep3Vector(0., tpToBase + tpDy/2, tpDz/2);
          _RightWallSheets.push_back(top);

          // Both are one plate thick, so either gives the step.
          rwAt += rwDepthDir*lwDx;

        } else if (kind == "Cu") {
          const double dx = _config.getDouble("stm.RightWall.copperSheet.dx");
          const double dy = _config.getDouble("stm.RightWall.copperSheet.dy");
          const double dz = _config.getDouble("stm.RightWall.copperSheet.dz");
          const double toBase =
            _config.getDouble("stm.RightWall.copperSheet.bottomToBase");
          // This one does not start at the reference: it begins a
          // lining's thickness downstream, clear of the front
          // shielding's copper.
          const double fromRefZ =
            _config.getDouble("stm.RightWall.copperSheet.fromRefZ");

          WallSheet s;
          s.name     = "CopperRightWallSheetPV";
          s.material = _config.getString("stm.RightWall.copperSheet.material");
          s.halfDim  = CLHEP::Hep3Vector(dx/2, dy/2, dz/2);
          s.center   = rwAt + rwDepthDir*(dx/2)
                     + CLHEP::Hep3Vector(0., toBase + dy/2, fromRefZ + dz/2);
          _RightWallSheets.push_back(s);

          rwAt += rwDepthDir*dx;

        } else {
          throw cet::exception("GEOM")
            << "STMMaker: " << key.str() << " is \"" << kind
            << "\", which is not one of Pb, BPouter, BPinner or Cu.\n";
        }
      }

    } else {
      _RightShieldingBuild  = _config.getBool("stm.RightShielding.build");
      _RightS_Length        = _config.getDouble("stm.RightShielding.Length");
      _RightSleaddepth      = _config.getDouble("stm.RightShielding.leaddepth");
      _RightScopperdepth    = _config.getDouble("stm.RightShielding.copperdepth");
      _RightSBPdepth        = _config.getDouble("stm.RightShielding.BPdepth");
      _RightSXmax           = _config.getDouble("stm.RightShielding.Right_Xmax");
    }

    // The updated top wall replaces TopShielding, as above.
    if (_handstacked) {
      // ---- Top wall ---------------------------------------------------
      //
      // The roof: three layers worked DOWNWARD in -y from the top of
      // the stack -- a pair of borated poly sheets, a layer of lead
      // bricks laid flat, and a pair of aluminium plates.
      //
      // The reference is inboard of the house reference in x, at that
      // reference's z, and at the TOP of the stack in y. The height is
      // built from the aluminium upward: wallToBase is the plates'
      // underside, and the reference sits the three layers' thickness
      // above it, so thickening a layer lifts the roof rather than
      // driving the aluminium into the wall below.

      _TopWallBuild = _config.getBool("stm.TopWall.build");

      const double twRefX =
        _ShieldHouseRefX + _config.getDouble("stm.TopWall.refFromRightWallX");
      const double twRefZ = _ShieldHouseRefZ;

      const std::string twPolyMat = _config.getString("stm.TopWall.BPmaterial");

      // The layers' own thicknesses, needed before any is placed so
      // that the reference can sit on top of the finished stack.
      const double twPolyT  = _config.getDouble("stm.TopWall.polySheet1.dy");
      const double twPlateT = _config.getDouble("stm.TopWall.plate.dy");

      // Courses run side by side across the roof and each course runs
      // along z. Both step away from the reference, so the directions
      // are negative and the layer hangs below it.
      const CLHEP::Hep3Vector twCourseDir( 0.,  0., 1.);
      const CLHEP::Hep3Vector twPitchDir (-1.,  0., 0.);
      const CLHEP::Hep3Vector twDepthDir ( 0., -1., 0.);

      // The bricks lie flat, so the layer's thickness is a brick's 2 in
      // whichever way the course is turned. Read from the first course,
      // and needed here rather than from layLayer because the roof's
      // height is built from the layer thicknesses before any brick is
      // placed.
      std::vector<int> twCourse1Types;
      _config.getVectorInt("stm.TopWall.leadLayer1Course1", twCourse1Types);
      if (twCourse1Types.empty()) {
        throw cet::exception("GEOM")
          << "STMMaker: stm.TopWall.leadLayer1Course1 is an empty course.\n";
      }
      const double twLeadT =
        brickSpan(_leadBrickDims, twCourse1Types[0],
                  rotationFor(_config.getString(
                    "stm.TopWall.leadLayer1Course1Orientation")),
                  twDepthDir, "stm.TopWall.leadLayer1Course1");

      const double twRefY = -_STM_SSCboreToBase
                          + _config.getDouble("stm.TopWall.wallToBase")
                          + twPlateT + twLeadT + twPolyT;

      // Walk the layers downward. Each starts at the face of the one
      // before it, so a layer's own thickness is the only thing that
      // separates it from the next. As in the other three walls the
      // cursor is a point advanced along depthDir, so the placements
      // below read the same way here as in the wall that stacks up.
      const int twLayerN = _config.getInt("stm.TopWall.layerN");
      CLHEP::Hep3Vector twAt(twRefX, twRefY, twRefZ);

      for (int L = 1; L <= twLayerN; ++L) {
        std::ostringstream key;
        key << "stm.TopWall.layer" << L;
        const std::string kind = _config.getString(key.str());

        if (kind == "BP") {
          // Two sheets side by side across the roof, the second also
          // starting a stated distance along z from the reference.
          double acrossX = 0.;
          for (int s = 1; s <= 2; ++s) {
            std::ostringstream base;
            base << "stm.TopWall.polySheet" << s;
            const double dx = _config.getDouble(base.str() + ".dx");
            const double dy = _config.getDouble(base.str() + ".dy");
            const double dz = _config.getDouble(base.str() + ".dz");
            // The first sheet starts at the reference; only the second
            // states an offset.
            const double fromRefZ = (s == 1) ? 0.
              : _config.getDouble(base.str() + ".fromRefZ");

            std::ostringstream name;
            name << "BPTopWallSheet" << s << "PV";

            WallSheet sheet;
            sheet.name     = name.str();
            sheet.material = twPolyMat;
            sheet.halfDim  = CLHEP::Hep3Vector(dx/2, dy/2, dz/2);
            sheet.center   = twAt
                           + twPitchDir*(acrossX + dx/2)
                           + twDepthDir*(dy/2)
                           + twCourseDir*(fromRefZ + dz/2);
            _TopWallSheets.push_back(sheet);

            acrossX += dx;
          }

          twAt += twDepthDir*twPolyT;

        } else if (kind == "Pb") {
          // Five courses side by side across the roof, each running
          // along z from its own start. Orientation and z offset are
          // per course, since the first course is turned differently
          // from the rest and the last two begin further along -- so no
          // layer-wide orientation is passed.
          //
          // The one lead layer, hence leadLayer1 rather than a counter.
          // The layer hangs below its top, which twDepthDir carries.
          const std::vector<BrickWallBrick> bricks =
            layLayer(_config, "stm.TopWall.leadLayer1", twAt,
                     twCourseDir, twPitchDir, twDepthDir,
                     _leadBrickDims, "");

          // The roof's height was built from twLeadT, read off the first
          // brick of course 1, before any brick was laid. Every brick in
          // the layer has to be that thick in y, or the layer would poke
          // into the plates below or the poly above.
          for (auto const & b : bricks) {
            const double t = brickSpan(_leadBrickDims, b.type,
                                       rotationFor(b.orientation), twDepthDir,
                                       "stm.TopWall.leadLayer1");
            if (std::abs(t - twLeadT) > 0.001) {
              throw cet::exception("GEOM")
                << "STMMaker: a stm.TopWall.leadLayer1 brick (type " << b.type
                << ", orientation " << b.orientation << ") is " << t
                << " mm thick in y, but the roof's height was built from "
                << twLeadT << " mm, the thickness of course 1's first brick.\n";
            }
          }

          _TopWallLeadLayers.push_back(
            BrickWall(_TopWallBuild, twAt,
                      twCourseDir, twPitchDir, twDepthDir,
                      bricks, std::vector<BrickWallBore>()));

          twAt += twDepthDir*twLeadT;

        } else if (kind == "Al") {
          // Two identical plates side by side in z.
          const std::string mat = _config.getString("stm.TopWall.plate.material");
          const double dx = _config.getDouble("stm.TopWall.plate.dx");
          const double dy = _config.getDouble("stm.TopWall.plate.dy");
          const double dz = _config.getDouble("stm.TopWall.plate.dz");
          const int    n  = _config.getInt("stm.TopWall.plateN");

          for (int p = 0; p < n; ++p) {
            std::ostringstream name;
            name << "AluminumTopWallPlate" << p+1 << "PV";

            WallSheet plate;
            plate.name     = name.str();
            plate.material = mat;
            plate.halfDim  = CLHEP::Hep3Vector(dx/2, dy/2, dz/2);
            plate.center   = twAt
                           + twPitchDir*(dx/2)
                           + twDepthDir*(dy/2)
                           + twCourseDir*(p*dz + dz/2);
            _TopWallSheets.push_back(plate);
          }

          twAt += twDepthDir*twPlateT;

        } else {
          throw cet::exception("GEOM")
            << "STMMaker: " << key.str() << " is \"" << kind
            << "\", which is not one of BP, Pb or Al.\n";
        }
      }

    } else {
      _TopShieldingBuild      = _config.getBool("stm.TopShielding.build");
      _TopShieldingSkirtBuild = _config.getBool("stm.TopShielding.Skirtbuild");
      _TopLiftBeam_L        = _config.getDouble("stm.TopShielding.LiftBeam_L");
      _TopLiftBeam_H        = _config.getDouble("stm.TopShielding.LiftBeam_H");
      _TopLiftBeam_T        = _config.getDouble("stm.TopShielding.LiftBeam_T");
      _TopLiftBeam_Xmove    = _config.getDouble("stm.TopShielding.LiftBeam_Xmove");
      _TopSZlength          = _config.getDouble("stm.TopShielding.Zlength");
      _TopSXlength          = _config.getDouble("stm.TopShielding.Xlength");
      _TopSFront_LT         = _config.getDouble("stm.TopShielding.Front_LT");
      _TopTFZlength         = _config.getDouble("stm.TopShielding.TFZlength");
      _TopTFXlength         = _config.getDouble("stm.TopShielding.TFXlength");
      _TopTBZlength         = _config.getDouble("stm.TopShielding.TBZlength");
      _TopScontainerdepth   = _config.getDouble("stm.TopShielding.containerdepth");
      _TopSleaddepth        = _config.getDouble("stm.TopShielding.leaddepth");
      _TopScopperdepth      = _config.getDouble("stm.TopShielding.copperdepth");
      _TopSBPdepth          = _config.getDouble("stm.TopShielding.BPdepth");
      _TopSZHole            = _config.getDouble("stm.TopShielding.Zhole");
      _TopSBarLeft          = _config.getDouble("stm.TopShielding.BarLeft");
      _TopSBarRight         = _config.getDouble("stm.TopShielding.BarRight");
      _TopSGapLeft          = _config.getDouble("stm.TopShielding.GapLeft");
      _TopSGapRight         = _config.getDouble("stm.TopShielding.GapRight");
      _TopSLeak             = _config.getDouble("stm.TopShielding.Leak");
    }

    _BackShieldingBuild   = _config.getBool("stm.BackShielding.build");
    _BackSBPThick         = _config.getDouble("stm.BackShielding.BPThick");
    _BackSBPLength        = _config.getDouble("stm.BackShielding.BPLength");
    _BackSBPHeight        = _config.getDouble("stm.BackShielding.BPHeight");
    _BackS_dX             = _config.getDouble("stm.BackShielding.BackS_dX");
    _BackS_dY             = _config.getDouble("stm.BackShielding.BackS_dY");
    _BackSPipeGap         = _config.getDouble("stm.BackShielding.ShieldingPipeGap");
    // Earlier geometries hard-code "BP", so only the hand-stacked one
    // reads the material.
    if (_handstacked) {
      _BackSBPMaterial    = readStringIfBuilt(_config, "stm.BackShielding.BPmaterial",
                                              _BackShieldingBuild);
    }

    // FrontToWall against the house as built.
    //
    // The earlier description placed the back poly off the house's own
    // length, so when FrontToWall disagreed with it the poly stopped short
    // of the hall wall, leaving space downstream of it (122.1 mm in
    // STM_v09). The updated one places the poly off FrontToWall instead
    // (see constructSTM), so its back face is always on the hall wall --
    // and a disagreement moves into the cable gap: lead to poly becomes
    // FrontToWall + ShieldingPipeGap - (where the lead actually ends)
    // rather than ShieldingPipeGap. This holds FrontToWall to the house.
    //
    // FrontToWall is defined against the back of the house LEAD, so only
    // the lead of the four walls is held to it; sheets and plates may
    // reach past it, as the top wall's aluminium plates do by design.
    if (_handstacked) {
      const double dimTolerance = 0.001;  // mm
      const CLHEP::Hep3Vector zDir(0., 0., 1.);

      double leadBackZ = -1.e30;
      for (auto const * layers : {&_RightWallLeadLayers, &_LeftWallLeadLayers,
                                  &_TopWallLeadLayers,   &_BottomWallLeadLayers}) {
        for (auto const & layer : *layers) {
          for (auto const & b : layer.bricks()) {
            const double half = brickSpan(_leadBrickDims, b.type,
                                          rotationFor(b.orientation), zDir,
                                          "the shield house") / 2.;
            leadBackZ = std::max(leadBackZ, b.center.z() + half);
          }
        }
      }

      const double gap = _STM_SSCFrontToWall + _BackSPipeGap - leadBackZ;

      if (leadBackZ > _STM_SSCFrontToWall + dimTolerance) {
        throw cet::exception("GEOM")
          << "STMMaker: the shield house lead ends at z = " << leadBackZ
          << " mm from the cradle front, "
          << leadBackZ - _STM_SSCFrontToWall << " mm past"
          << " stm.STM_SSC.FrontToWall (" << _STM_SSCFrontToWall << " mm).\n"
          << "The cable gap behind the lead is " << gap
          << " mm instead of ShieldingPipeGap (" << _BackSPipeGap << " mm)."
          << " Check FrontToWall against the walls' reference z and course"
          << " lengths.\n";
      }
      else if (leadBackZ < _STM_SSCFrontToWall - dimTolerance) {
        mf::LogWarning("GEOM")
          << "STMMaker: the shield house lead ends at z = " << leadBackZ
          << " mm, " << _STM_SSCFrontToWall - leadBackZ << " mm short of"
          << " stm.STM_SSC.FrontToWall (" << _STM_SSCFrontToWall << " mm), so"
          << " the cable gap behind it is " << gap
          << " mm instead of ShieldingPipeGap (" << _BackSPipeGap << " mm).\n";
      }
      else if (_verbosityLevel > 0) {
        mf::LogInfo("GEOM")
          << "STMMaker: the shield house lead ends at FrontToWall (<"
          << dimTolerance << " mm), at z = " << leadBackZ << " mm.\n";
      }
    }


    _InnerShieldingBuild        = _config.getBool("stm.InnerShielding.build");

    // The updated inner shielding. Unlike the walls this is not a
    // stack of layers: every piece is placed on its own, so what
    // follows is three lists rather than a sequence.
    //
    // The reference is the bottom wall's, moved in 0.5 in and up 6 in
    // -- directly above the origin of that wall's first copper prism,
    // and level with the top of its lead.
    if (_handstacked) {
      const double innerRefX = _ShieldHouseRefX
                             + _config.getDouble("stm.BottomWall.refFromHouseX")
                             + _config.getDouble("stm.InnerShielding.refFromBottomWallX");
      const double innerRefY = -_STM_SSCboreToBase
                             + _config.getDouble("stm.InnerShielding.refFromBottomWallY");
      const double innerRefZ = _ShieldHouseRefZ;

      const std::string innerCuMat =
        _config.getString("stm.InnerShielding.lining.material");

      // The holes through this section. Stated here rather than taken
      // from the front shielding so that the section stands on its own:
      // the bore follows the collimator either way, and this one would
      // otherwise be unbuildable with the front shielding switched off.
      _InnerShieldingBores =
        readBores(_config, "stm.InnerShielding", _STM_SSCoffset_Spot);
      const size_t innerBoreN = _InnerShieldingBores.size();

      // Which bores a named piece is drilled for.
      auto innerBores = [&](std::string const & key) {
        return readPieceBores(_config,
                              "stm.InnerShielding." + key + ".bores",
                              innerBoreN,
                              "stm.InnerShielding." + key);
      };

      // ---- the copper lining ----------------------------------------
      //
      // The left slanted piece starts at the reference; the other five
      // start a step above it. Every piece here is pinned at the
      // mid-plane of its own height, since an extruded solid is centered
      // on its placement point and the boxes follow the same height.
      const double innerStepUp   = _config.getDouble("stm.InnerShielding.lining.stepUp");
      const double innerCuHeight = _config.getDouble("stm.InnerShielding.lining.height");

      // The mid-plane the five stepped-up pieces share.
      const double innerCuY = innerRefY + innerStepUp + innerCuHeight/2;

      // One helper for the four lining prisms and the four lead ones:
      // they differ only in which keys they read.
      auto innerPrism = [&](std::string const & key,
                            std::string const & name,
                            std::string const & material,
                            double length,
                            double y) {
        InnerShieldingPrism p;
        p.name     = name;
        p.material = material;
        _config.getVectorDouble("stm.InnerShielding." + key + ".UVerts", p.uVerts);
        _config.getVectorDouble("stm.InnerShielding." + key + ".VVerts", p.vVerts);
        p.length      = length;
        p.orientation = _config.getString("stm.InnerShielding." + key + ".orientation");
        p.anchor      = CLHEP::Hep3Vector(
          innerRefX + _config.getDouble("stm.InnerShielding." + key + ".fromRefX"),
          y,
          innerRefZ + _config.getDouble("stm.InnerShielding." + key + ".fromRefZ"));
        // The key says which holes it takes, so the caller no longer
        // has to; a piece that names none takes none.
        p.bores = innerBores(key);
        return p;
      };

      // The left slanted piece: the one at the reference height.
      {
        const double h = _config.getDouble("stm.InnerShielding.liningLeftSlanted.height");
        _InnerShieldingPrisms.push_back(
          innerPrism("liningLeftSlanted", "CopperInnerLeftSlantedPV", innerCuMat,
                     h, innerRefY + h/2));
      }

      // The other three lining prisms, a step up from the reference.
      _InnerShieldingPrisms.push_back(
        innerPrism("liningMidLeftWall", "CopperInnerMidLeftWallPV", innerCuMat,
                   innerCuHeight, innerCuY));
      _InnerShieldingPrisms.push_back(
        innerPrism("liningMidSlanted", "CopperInnerMidSlantedPV", innerCuMat,
                   innerCuHeight, innerCuY));
      _InnerShieldingPrisms.push_back(
        innerPrism("liningMidBack", "CopperInnerMidBackPV", innerCuMat,
                   innerCuHeight, innerCuY));

      // The two lining boxes, on that same mid-plane.
      {
        const double rdx = _config.getDouble("stm.InnerShielding.liningMidRight.dx");
        const double rdz = _config.getDouble("stm.InnerShielding.liningMidRight.dz");
        InnerShieldingBox right;
        right.name        = "CopperInnerMidRightPV";
        right.material    = innerCuMat;
        right.halfDim     = CLHEP::Hep3Vector(rdx/2, innerCuHeight/2, rdz/2);
        right.center      = CLHEP::Hep3Vector(
          innerRefX + _config.getDouble("stm.InnerShielding.liningMidRight.centerFromRefX"),
          innerCuY,
          innerRefZ + _config.getDouble("stm.InnerShielding.liningMidRight.centerFromRefZ"));
        right.orientation = "000";
        right.bores       = innerBores("liningMidRight");
        right.nudge       = false;
        _InnerShieldingBoxes.push_back(right);

        const double bdx = _config.getDouble("stm.InnerShielding.liningMidBackBox.dx");
        const double bdz = _config.getDouble("stm.InnerShielding.liningMidBackBox.dz");
        InnerShieldingBox back;
        back.name        = "CopperInnerMidBackBoxPV";
        back.material    = innerCuMat;
        back.halfDim     = CLHEP::Hep3Vector(bdx/2, innerCuHeight/2, bdz/2);
        back.center      = CLHEP::Hep3Vector(
          innerRefX + _config.getDouble("stm.InnerShielding.liningMidBackBox.centerFromRefX"),
          innerCuY,
          innerRefZ + _config.getDouble("stm.InnerShielding.liningMidBackBox.centerFromRefZ"));
        back.orientation = "000";
        back.bores       = innerBores("liningMidBackBox");
        back.nudge       = false;
        _InnerShieldingBoxes.push_back(back);
      }

      // ---- the lead -------------------------------------------------

      // The bricks that lie square, each given by its own center as
      // {type, orientation, x, y, z}.
      {
        const int n = _config.getInt("stm.InnerShielding.leadBrickN");
        for (int i = 1; i <= n; ++i) {
          std::ostringstream key;
          key << "stm.InnerShielding.leadBrick" << i;
          std::vector<double> v;
          _config.getVectorDouble(key.str(), v, 5);

          const int type = int(v[0]);
          if (type < 1 || type > int(_leadBrickDims.size())) {
            throw cet::exception("GEOM")
              << "STMMaker: " << key.str() << " names brick type " << type
              << ", but stm.leadBrick.typeN defines only "
              << _leadBrickDims.size() << ".\n";
          }
          // The orientation rides in the tuple as a number, so it is
          // turned back into the three-digit code the resolver reads;
          // leading zeros matter, hence the width.
          std::ostringstream orient;
          orient << std::setw(3) << std::setfill('0') << int(v[1]);

          BrickWallBrick b;
          b.type        = type;
          b.orientation = orient.str();
          b.center      = CLHEP::Hep3Vector(innerRefX + v[2],
                                            innerRefY + v[3],
                                            innerRefZ + v[4]);
          // The one on the LaBr axis names that bore; the rest state
          // no bores key and take none.
          std::ostringstream bkey;
          bkey << "leadBrick" << i;
          b.bores       = innerBores(bkey.str());
          _InnerShieldingBricks.push_back(b);
        }
      }

      // The three turned 45 degrees about y. Stated as boxes because
      // the 45 degree code is a special case in OrientationResolver
      // that cannot follow a 90 degree step.
      //
      // Their nominal sizes would touch the volumes around them once
      // turned, so each half-dimension is backed off by the shared
      // nudge rather than the config carrying shrunk numbers.
      {
        const std::string orient =
          _config.getString("stm.InnerShielding.leadAngledBox.orientation");
        for (int i = 1; i <= 3; ++i) {
          std::ostringstream base, name;
          base << "stm.InnerShielding.leadAngledBox" << i;
          name << "LeadInnerAngledBox" << i << "PV";

          const double dx = _config.getDouble(base.str() + ".dx");
          const double dy = _config.getDouble(base.str() + ".dy");
          const double dz = _config.getDouble(base.str() + ".dz");

          InnerShieldingBox b;
          b.name        = name.str();
          b.material    = _leadBrickMaterial;
          // Nominal. The construction code backs each half-dimension
          // off by its own nudge, since the constant lives there.
          b.halfDim     = CLHEP::Hep3Vector(dx/2, dy/2, dz/2);
          b.center      = CLHEP::Hep3Vector(
            innerRefX + _config.getDouble(base.str() + ".centerFromRefX"),
            innerRefY + _config.getDouble(base.str() + ".centerFromRefY"),
            innerRefZ + _config.getDouble(base.str() + ".centerFromRefZ"));
          b.orientation = orient;
          std::ostringstream bkey;
          bkey << "leadAngledBox" << i;
          b.bores       = innerBores(bkey.str());
          b.nudge       = true;
          _InnerShieldingBoxes.push_back(b);
        }
      }

      // The angled triangular prism.
      _InnerShieldingPrisms.push_back(
        innerPrism("leadAngled", "LeadInnerAngledPV", _leadBrickMaterial,
                   _config.getDouble("stm.InnerShielding.leadAngled.length"),
                   innerRefY + _config.getDouble("stm.InnerShielding.leadAngled.fromRefY")));

      // The bored trapezoid.
      _InnerShieldingPrisms.push_back(
        innerPrism("leadHoleCut", "LeadInnerHoleCutPV", _leadBrickMaterial,
                   _config.getDouble("stm.InnerShielding.leadHoleCut.length"),
                   innerRefY + _config.getDouble("stm.InnerShielding.leadHoleCut.fromRefY")));

      // The two clipped prisms, above and below the trapezoid. They
      // share one outline and differ only in thickness and height, so
      // the shared keys are read once.
      {
        std::vector<double> u, v;
        _config.getVectorDouble("stm.InnerShielding.leadClipped.UVerts", u);
        _config.getVectorDouble("stm.InnerShielding.leadClipped.VVerts", v);
        const std::string orient =
          _config.getString("stm.InnerShielding.leadClipped.orientation");

        for (int i = 1; i <= 2; ++i) {
          std::ostringstream base, name;
          base << "stm.InnerShielding.leadClipped" << i;
          name << "LeadInnerClipped" << i << "PV";

          InnerShieldingPrism p;
          p.name        = name.str();
          p.material    = _leadBrickMaterial;
          p.uVerts      = u;
          p.vVerts      = v;
          p.length      = _config.getDouble(base.str() + ".length");
          p.orientation = orient;
          p.anchor      = CLHEP::Hep3Vector(
            innerRefX + _config.getDouble(base.str() + ".fromRefX"),
            innerRefY + _config.getDouble(base.str() + ".fromRefY"),
            innerRefZ + _config.getDouble(base.str() + ".fromRefZ"));
          std::ostringstream bkey;
          bkey << "leadClipped" << i;
          p.bores       = innerBores(bkey.str());
          _InnerShieldingPrisms.push_back(p);
        }
      }

      // ---- bore fit check ------------------------------------------
      //
      // As for the front shielding's bricks: every hole a piece names
      // must lie inside that piece across the beam, or the subtraction
      // clips only part of the hole -- or nothing -- and the beam sees
      // solid material with nothing to say so. Each piece's x and y
      // extent below is the one the hole has to clear over the piece's
      // whole depth along the beam.
      const CLHEP::Hep3Vector xDir(1., 0., 0.), yDir(0., 1., 0.), zDir(0., 0., 1.);
      auto checkFit = [&](std::string const & what,
                          double xlo, double xhi, double ylo, double yhi,
                          std::vector<int> const & bores) {
        for (int id : bores) {
          BrickWallBore const & b = _InnerShieldingBores.at(id-1);
          if (b.center.x() - b.radius < xlo || b.center.x() + b.radius > xhi ||
              b.center.y() - b.radius < ylo || b.center.y() + b.radius > yhi) {
            throw cet::exception("GEOM")
              << "STMMaker: bore " << id << " (" << b.axis << ", radius "
              << b.radius << " at x = " << b.center.x() << ", y = "
              << b.center.y() << ") does not fit " << what
              << ", which spans x [" << xlo << ", " << xhi << "] and y ["
              << ylo << ", " << yhi << "] mm.\n";
          }
        }
      };

      // Bricks: the worn solid the hole is cut from. Each of the brick's
      // own axes lands on rot^-1 * e_k (as in spanAlong) and loses the
      // wear for whichever Mu2e direction it lands on.
      for (auto const & b : _InnerShieldingBricks) {
        if (b.bores.empty()) continue;
        const CLHEP::HepRotation rot = rotationFor(b.orientation);
        const CLHEP::Hep3Vector & d = _leadBrickDims.at(b.type-1);
        CLHEP::Hep3Vector worn;
        for (int k = 0; k < 3; ++k) {
          const CLHEP::Hep3Vector axis =
            rot.inverse() * CLHEP::Hep3Vector(k == 0, k == 1, k == 2);
          worn[k] = d[k] - 2.*((std::abs(axis.y()) > 0.5) ? _leadBrickWearY : _leadBrickWearXZ);
        }
        const double hx = spanAlong(worn, rot, xDir)/2.;
        const double hy = spanAlong(worn, rot, yDir)/2.;
        checkFit("a stm.InnerShielding lead brick",
                 b.center.x() - hx, b.center.x() + hx,
                 b.center.y() - hy, b.center.y() + hy, b.bores);
      }

      // Boxes: their nominal size, turned by their orientation.
      for (auto const & bx : _InnerShieldingBoxes) {
        if (bx.bores.empty()) continue;
        const CLHEP::HepRotation rot = rotationFor(bx.orientation);
        const double hx = spanAlong(2.*bx.halfDim, rot, xDir)/2.;
        const double hy = spanAlong(2.*bx.halfDim, rot, yDir)/2.;
        checkFit(bx.name, bx.center.x() - hx, bx.center.x() + hx,
                 bx.center.y() - hy, bx.center.y() + hy, bx.bores);
      }

      // Prisms: only those whose outline's v lands along the beam are
      // supported, which is every bored one here ("100"). The beam then
      // runs through the outline along v, so the hole must lie inside the
      // outline's u range at every v the piece covers -- the overlap of
      // those ranges -- and inside the sweep in the other direction.
      for (auto const & p : _InnerShieldingPrisms) {
        if (p.bores.empty()) continue;
        const CLHEP::HepRotation rotInv = rotationFor(p.orientation).inverse();
        const CLHEP::Hep3Vector uTo = rotInv * xDir;   // the outline's u
        const CLHEP::Hep3Vector vTo = rotInv * yDir;   // the outline's v
        const CLHEP::Hep3Vector sTo = rotInv * zDir;   // the sweep
        if (std::abs(vTo.z()) < 0.5 || std::abs(uTo.x()) < 0.5 ||
            std::abs(sTo.y()) < 0.5) {
          throw cet::exception("GEOM")
            << "STMMaker: " << p.name << " is bored, but its orientation "
            << p.orientation << " does not put its outline's v along the"
            << " beam, u along x and its sweep along y, which is all the"
            << " bore fit check handles.\n";
        }

        // The outline's u range on the line at a given v.
        const size_t n = p.uVerts.size();
        auto uRangeAt = [&](double v, double & lo, double & hi) {
          lo = 1.e30; hi = -1.e30;
          for (size_t i = 0; i < n; ++i) {
            const double u1 = p.uVerts[i],       v1 = p.vVerts[i];
            const double u2 = p.uVerts[(i+1)%n], v2 = p.vVerts[(i+1)%n];
            if ((v - v1)*(v - v2) > 0.) continue;    // edge misses this v
            if (v1 == v2) {                          // edge lies along it
              lo = std::min(lo, std::min(u1, u2));
              hi = std::max(hi, std::max(u1, u2));
            } else {                                 // edge crosses it
              const double u = u1 + (u2 - u1)*(v - v1)/(v2 - v1);
              lo = std::min(lo, u);
              hi = std::max(hi, u);
            }
          }
        };
        // The edges are straight, so the overlap is set at the vertices.
        double uLo = -1.e30, uHi = 1.e30;
        for (size_t i = 0; i < n; ++i) {
          double lo, hi;
          uRangeAt(p.vVerts[i], lo, hi);
          uLo = std::max(uLo, lo);
          uHi = std::min(uHi, hi);
        }

        const double sx = (uTo.x() > 0.) ? 1. : -1.;
        const double xlo = p.anchor.x() + std::min(sx*uLo, sx*uHi);
        const double xhi = p.anchor.x() + std::max(sx*uLo, sx*uHi);
        checkFit(p.name, xlo, xhi,
                 p.anchor.y() - p.length/2., p.anchor.y() + p.length/2.,
                 p.bores);
      }
    }



    _ElectronicShieldingBuild      = _config.getBool("stm.ElectronicShielding.build");
    _ElectronicSSiGridX    = _config.getDouble("stm.ElectronicShielding.SiGridX");
    _ElectronicSSiGridY    = _config.getDouble("stm.ElectronicShielding.SiGridY");
    _ElectronicSSiGridZ    = _config.getDouble("stm.ElectronicShielding.SiGridZ");
    if (!_handstacked) {
      // The earlier description: one block and one grid, placed from the
      // front shielding.
      _ElectronicSSiXcenter  = _config.getDouble("stm.ElectronicShielding.SiXcenter");
      _ElectronicSSiYcenter  = _config.getDouble("stm.ElectronicShielding.SiYcenter");
      _ElectronicSSiZcenter  = _config.getDouble("stm.ElectronicShielding.SiZcenter");
      _ElectronicSConcreteT  = _config.getDouble("stm.ElectronicShielding.ConcreteT");
      _ElectronicSGapToSi    = _config.getDouble("stm.ElectronicShielding.GapToSi");
    }
    else if (_ElectronicShieldingBuild) {
      // The hand-stacked description: lab standard concrete blocks stacked
      // on the hall floor and two silicon tile grids, all placed in the
      // hall rather than from the shield house. Resolved in make(), where
      // the hall's east wall is known.
      _ElectronicSXOffset             = _config.getDouble("stm.ElectronicShielding.xOffset");
      _ElectronicSConcreteDims        = CLHEP::Hep3Vector(
        _config.getDouble("stm.ElectronicShielding.concrete.dx"),
        _config.getDouble("stm.ElectronicShielding.concrete.dy"),
        _config.getDouble("stm.ElectronicShielding.concrete.dz"));
      _ElectronicSConcreteBlockN      = _config.getInt(   "stm.ElectronicShielding.concrete.blockN");
      _ElectronicSConcreteBlockGap    = _config.getDouble("stm.ElectronicShielding.concrete.blockGap");
      _ElectronicSConcreteFromGround  = _config.getDouble("stm.ElectronicShielding.concrete.fromGround");
      _ElectronicSConcreteToWall      = _config.getDouble("stm.ElectronicShielding.concrete.downstreamFaceToWall");
      _ElectronicSSiNX                = _config.getInt(   "stm.ElectronicShielding.si.nX");
      _ElectronicSSiNY                = _config.getInt(   "stm.ElectronicShielding.si.nY");
      _ElectronicSSiFromGround        = _config.getDouble("stm.ElectronicShielding.si.fromGround");
      _config.getVectorDouble("stm.ElectronicShielding.si.gapsToConcrete",
                              _ElectronicSSiGapsToConcrete);

      if (_ElectronicSConcreteBlockN < 1 || _ElectronicSSiNX < 1 || _ElectronicSSiNY < 1) {
        throw cet::exception("GEOM")
          << "STMMaker: stm.ElectronicShielding needs at least one concrete"
          << " block and one silicon tile in x and y (blockN = "
          << _ElectronicSConcreteBlockN << ", si.nX = " << _ElectronicSSiNX
          << ", si.nY = " << _ElectronicSSiNY << ").\n";
      }
    }


    // The SSC front shield replaces the absorber in the updated geometry,
    // so only one description's keys are read.
    if (!_handstacked) {
      _STM_AbsorberBuild          = _config.getBool("stm.STM_Absorber.build");
      _STM_Absorber_hW            = _config.getDouble("stm.STM_Absorber.hW");
      _STM_Absorber_hH            = _config.getDouble("stm.STM_Absorber.hH");
      _STM_Absorber_hT            = _config.getDouble("stm.STM_Absorber.hT");
      _STM_Absorber_GaptoSSC      = _config.getDouble("stm.STM_Absorber.GaptoSSC");
    }

  }
} // namespace mu2e
