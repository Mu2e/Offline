#ifndef ExtinctionMonitorFNAL_Utilities_inc_EMFBoxIO_hh
#define ExtinctionMonitorFNAL_Utilities_inc_EMFBoxIO_hh

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

    } // IO
  } // ExtMonFNAL
} // mu2e

#endif/*ExtinctionMonitorFNAL_Utilities_inc_EMFBoxIO_hh*/
