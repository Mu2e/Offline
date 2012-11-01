#ifndef ExtinctionMonitorFNAL_Utilities_inc_EMFBoxIO_hh
#define ExtinctionMonitorFNAL_Utilities_inc_EMFBoxIO_hh

#include "MCDataProducts/inc/MARSInfo.hh"

namespace mu2e {
  namespace ExtMonFNAL {
    namespace IO {

      struct EMFBoxHit {
        double emx;
        double emy;
        double emz;
        double mu2epx;
        double mu2epy;
        double mu2epz;
        double time;
        int pdgId;
        int vdId;

        EMFBoxHit() : emx(), emy(), emz(), mu2epx(), mu2epy(), mu2epz(), time(), pdgId(), vdId() {}

        static const char *branchDescription() {
          return "emx/D:emy/D:emz/D:mu2epx/D:mu2epy/D:mu2epz/D:time/D:pdgId/I:vdId/I";
        }
      };

      //================================================================
      struct EMFRoomHit {
        double dumpx;
        double dumpy;
        double dumpz;
        double mu2epx;
        double mu2epy;
        double mu2epz;
        double time;
        int pdgId;
        int srcType;

        EMFRoomHit() : dumpx(), dumpy(), dumpz(), mu2epx(), mu2epy(), mu2epz(), time(), pdgId(), srcType() {}

        static const char *branchDescription() {
          return "dumpx/D:dumpy/D:dumpz/D:mu2epx/D:mu2epy/D:mu2epz/D:time/D:pdgId/I:srcType/I";
        }
      };

      //================================================================
      struct ParticleRandomization {

        double sigmax;
        double sigmay;
        double sigmaz;
        double correlationCoefficient;  // between the two coordinates relevant for the source plane

        // "Gaussian" distribution of neighbor directions (nx', ny')
        // w.r.t this direction leads to Rayleigh distribution of
        // r = sin(theta'):  pdf ~ r/sigma^2 exp(-r^2/(2*sigma^2))
        //
        // ML estimate of sigma is sqrt(\sum_1^N x^2_i /(2N))

        double rSigmaML;

        ParticleRandomization()
          : sigmax(), sigmay(), sigmaz(), correlationCoefficient()
          , rSigmaML()
        {}

        static const char *branchDescription() {
          return "sigmax/D:sigmay/D:sigmaz/D:correlationCoefficient/D:rSigmaML/D";
        }
      };

      //================================================================
      struct StoppedMuon {
        double emx;
        double emy;
        double emz;
        double time;
        double endek;

        int pdgId;

        unsigned endG4Status;
        unsigned stoppingCode;

        int stoppedInSensor; // would bool be OK for ROOT?

        StoppedMuon() : emx(), emy(), emz(), time(), endek(), pdgId(), endG4Status(), stoppingCode(), stoppedInSensor() {}

        static const char *branchDescription() {
          return "emx/D:emy/D:emz/D:time/D:endek/D:pdgId/I:endG4Status/i:stoppingCode/i:stoppedInSensor/I";
        }
      };

      //================================================================
      struct MARSInfo {
        mu2e::MARSInfo info;
        static const char *branchDescription() {
          return "weight/D:protonNumber/I:subRunNumber/I:runNumber/I";
        }
      };

      //================================================================
      struct G4JobInfo {
        unsigned run;
        unsigned subrun;
        unsigned event;
        G4JobInfo() : run(), subrun(), event() {}
        G4JobInfo(unsigned r, unsigned s, unsigned e) : run(r), subrun(s), event(e) {}

        static const char *branchDescription() {
          return "run/i:subrun/i:event/i";
        }
      };

      // a "less than" comparison
      struct CmpG4JobInfo {
        bool operator()(const G4JobInfo& a, const G4JobInfo& b) const {
          return
            (a.run < b.run) || ((a.run == b.run) &&
                                ((a.subrun < b.subrun) || ((a.subrun == b.subrun) &&
                                                           (a.event < b.event))));
        }
      };

      bool operator==(const G4JobInfo& a, const G4JobInfo& b) {
        return (a.event == b.event) && (a.subrun == b.subrun) && (a.run == b.run);
      }

      //================================================================

    } // IO
  } // ExtMonFNAL
} // mu2e

#endif/*ExtinctionMonitorFNAL_Utilities_inc_EMFBoxIO_hh*/
