// Andrei Gaponenko, 2012

#include "ExtinctionMonitorFNAL/Utilities/inc/EMFRandomizationSourceDefs.hh"

#include <cmath>
#include <cstdlib>
#include <cassert>

#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"

namespace mu2e {
  namespace ExtMonFNAL {
    namespace Randomization {

      //================================================================
      SourcePlaneGeom::SourcePlaneGeom(const fhicl::ParameterSet& pset)
        : zFront(pset.get<double>("zFront"))
        , xSW(pset.get<double>("xSW"))
        , yFloor(pset.get<double>("yFloor"))
        , zBack(pset.get<double>("zBack"))
        , xNE(pset.get<double>("xNE"))
        , yCeiling(pset.get<double>("yCeiling"))
        , srcPositionTolerance(pset.get<double>("srcPositionTolerance"))
        , signalHalfdx(pset.get<double>("signalHalfdx"))
        , signalHalfdy(pset.get<double>("signalHalfdy"))
      {
        if(xSW >= xNE) {
          throw cet::exception("BADCONFIG")
            <<"Error: SourcePlaneGeom: xSW ("<<xSW<<") >= xNE ("<<xNE<<")";
        }

        if(yFloor >= yCeiling) {
          throw cet::exception("BADCONFIG")
            <<"Error: SourcePlaneGeom: yFloor ("<<yFloor<<") >= yCeiling ("<<yCeiling<<")";
        }

        if(zBack >= zFront) {
          throw cet::exception("BADCONFIG")
            <<"Error: SourcePlaneGeom: zBack ("<<zBack<<") >= zFront ("<<zFront<<")";
        }
      }

      //================================================================
      std::ostream& operator<<(std::ostream&os, const SourcePlaneGeom& gm) {
        return
          os<<"SourcePlaneGeom(zFront="<<gm.zFront
            <<", xSW="<<gm.xSW
            <<", zBack="<<gm.zBack
            <<", yFloor="<<gm.yFloor
            <<", xNE="<<gm.xNE
            <<", yCeiling="<<gm.yCeiling
            <<", tolerance="<<gm.srcPositionTolerance
            <<", signalHalfdx="<<gm.signalHalfdx
            <<", signalHalfdy="<<gm.signalHalfdy
            <<" )"
          ;
      }

      //================================================================
      SourceType classifySource(const CLHEP::Hep3Vector& posMu2e,
                                const SourcePlaneGeom& srcGeom,
                                const ProtonBeamDump& dump,
                                const ExtMon& extmon)
      {
        const CLHEP::Hep3Vector posDump = dump.mu2eToBeamDump_position(posMu2e);

        // The order of the tests affects the result only for corner
        // cases, where classification is ambiguous.

        if(std::abs(posDump.z() - srcGeom.zFront) < srcGeom.srcPositionTolerance) {
          // check for Signal
          if(isSignal(extmon.mu2eToExtMon_position(posMu2e), srcGeom)) {
            return SourceSignal;
          }
          return SourceFront;
        }

        if(std::abs(posDump.x() - srcGeom.xSW) < srcGeom.srcPositionTolerance) return SourceSouthWest;
        if(std::abs(posDump.y() - srcGeom.yFloor) < srcGeom.srcPositionTolerance) return SourceFloor;
        if(std::abs(posDump.x() - srcGeom.xNE) < srcGeom.srcPositionTolerance) return SourceNorthEast;
        if(std::abs(posDump.y() - srcGeom.yCeiling) < srcGeom.srcPositionTolerance) return SourceCeiling;

        if(std::abs(posDump.z() - srcGeom.zBack) < srcGeom.srcPositionTolerance) return SourceBack;

        throw cet::exception("BADINPUTS")
          <<"Error: failed to assign input posDump = "<<posDump<<" to an input source plane\n";
        //return NUM_SOURCES;
      }

      //================================================================
      bool isSignal(const CLHEP::Hep3Vector& posExtMon,
                    const SourcePlaneGeom& srcGeom)
      {
        return
          (std::abs(posExtMon.x()) < srcGeom.signalHalfdx) &&
          (std::abs(posExtMon.y()) < srcGeom.signalHalfdy);
      }

      //================================================================
      bool inRange(SourceType st,
                   const CLHEP::Hep3Vector& posDump,
                   const SourcePlaneGeom& srcGeom,
                   const ProtonBeamDump& dump,
                   const ExtMon& extmon)
      {
        switch(st) {

        default: assert(false);

        case SourceFront: case SourceBack:
          return
            (srcGeom.xSW <= posDump.x()) && (posDump.x() <= srcGeom.xNE) &&
            (srcGeom.yFloor <= posDump.y()) && (posDump.y() <= srcGeom.yCeiling);

        case SourceSouthWest: case SourceNorthEast:
          return
            (srcGeom.zBack <= posDump.z()) && (posDump.z() <= srcGeom.zFront) &&
            (srcGeom.yFloor <= posDump.y()) && (posDump.y() <= srcGeom.yCeiling);

        case SourceFloor: case SourceCeiling:
          return
            (srcGeom.xSW <= posDump.x()) && (posDump.x() <= srcGeom.xNE) &&
            (srcGeom.zBack <= posDump.z()) && (posDump.z() <= srcGeom.zFront);

        case SourceSignal: {
          const CLHEP::Hep3Vector posExtMon =
            extmon.mu2eToExtMon_position(dump.beamDumpToMu2e_position(posDump));

          return isSignal(posExtMon, srcGeom);
        }
        } // switch(st)

      } // inRange()

      //================================================================
    } // namespace Randomization
  } // namespace ExtMonFNAL
} // namespace mu2e
