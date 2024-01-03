// Add random noise hits with the specified probability.
// Time stamps of created hits are in [0, numClockTicksPerDebuncherPeriod-1].
//
// Andrei Gaponenko, 2012

#ifndef ExtinctionMonitorFNAL_Digitization_PixelNoise_hh
#define ExtinctionMonitorFNAL_Digitization_PixelNoise_hh

#include "Offline/RecoDataProducts/inc/ExtMonFNALRawHit.hh"
#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/RandFlat.h"

#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

namespace mu2e {

  namespace ExtMonFNAL {

    class ExtMon;

    class PixelNoise {
    public:

      PixelNoise(art::RandomNumberGenerator::base_engine_t& rng,
                 const ExtMon **em,
                 const int numBCs,
                 double noisePerPixelPerBC);

      // Adds noise hits to the collection.  Noise hits can be
      // identified by the lack of associated truth.  No need to have a
      // truth collection ptr arg.
      void add(ExtMonFNALRawHitCollection *hits);

    private:
      CLHEP::RandPoissonQ poisson_;
      CLHEP::RandFlat flat_;
      const ExtMon **extmon_;
      int numBCs_;
      double noisePerPixelPerBC_;
    };

  } // namespace ExtMonFNAL
} // namespace mu2e

#endif/*ExtinctionMonitorFNAL_Digitization_PixelNoise_hh*/
