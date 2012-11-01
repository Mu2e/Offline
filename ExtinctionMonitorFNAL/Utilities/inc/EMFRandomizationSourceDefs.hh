// Andrei Gaponenko, 2012

#ifndef ExtinctionMonitorFNAL_Utilities_inc_EMFMARSRoomSourceDefs_hh
#define ExtinctionMonitorFNAL_Utilities_inc_EMFMARSRoomSourceDefs_hh

#include <ostream>

namespace fhicl { class ParameterSet; }
namespace CLHEP { class Hep3Vector; }

namespace mu2e {

  class ProtonBeamDump;

  namespace ExtMonFNAL {

    class ExtMon;

    namespace Randomization {

      //================================================================
      // src % 3 == 0: don't smear dumpz
      // src % 3 == 1: don't smear dumpx
      // src % 3 == 2: don't smear dumpy
      enum SourceType {
        SourceFront = 0,
        SourceSouthWest,
        SourceFloor,
        SourceBack,
        SourceNorthEast,
        SourceCeiling,

        SourceSignal, // SourceFront particles which are also signal candidates

        NUM_SOURCES
      };

      //================================================================
      // These are in beam dump coordinates
      // The numbers are used to
      //
      //    1) assign input MARS particles to a source plane
      //    2) establish boundaries which (non-signal) smearing should not cross

      struct SourcePlaneGeom {
        double zFront;
        double xSW;
        double yFloor;
        double zBack;
        double xNE;
        double yCeiling;

        double srcPositionTolerance;

        // Signal box in ExtMon coordinates is [-dx, +dx]*[-dy,+dy]
        double signalHalfdx;
        double signalHalfdy;

        SourcePlaneGeom(const fhicl::ParameterSet& pset);
      };


      std::ostream& operator<<(std::ostream&os, const SourcePlaneGeom& gm);

      //================================================================
      SourceType classifySource(const CLHEP::Hep3Vector& posMu2e,
                                const SourcePlaneGeom& srcGeom,
                                const ProtonBeamDump& dump,
                                const ExtMon& extmon);


      bool isSignal(const CLHEP::Hep3Vector& posExtMon, const SourcePlaneGeom& srcGeom);

      // This one is used to make sure a particle on a source plane
      // stays there after randomization - only the 2D in-plane position is checked.
      bool inRange(SourceType st,
                   const CLHEP::Hep3Vector& posDump,
                   const SourcePlaneGeom& srcGeom,
                   const ProtonBeamDump& dump,
                   const ExtMon& extmon);

      //================================================================
    }
  }
}

#endif/*ExtinctionMonitorFNAL_Utilities_inc_EMFMARSRoomSourceDefs_hh*/
